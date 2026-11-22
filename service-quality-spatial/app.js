'use strict';
/* ==========================================================================
 * 片区服务质量空间分析
 *  - 虚构 8×8 网格 + 虚构服务记录（固定生成种子，数据稳定）
 *  - 每片区：样本量 / 均值 / 有效覆盖
 *  - 邻接：共享边(rook) / 共享角(queen)
 *  - 全局 Moran's I + 局部 LISA，可重复随机置换（支持同类型片区内置换）
 *  - 局部检验经 Benjamini–Hochberg 多重比较校正
 *  - 无邻居 / 无记录 / 覆盖不足片区单列，绝不当作零值
 * ========================================================================== */

/* ---------------- 常量 ---------------- */
const GRID = 8;            // 8×8 网格
const DAYS = 14;           // 观察窗口：14 天
const MIN_N = 30;          // 参与计算的最小样本量
const MIN_COV = 0.5;       // 参与计算的最小有效覆盖（有天记录的天数占比）
const ALPHA = 0.05;        // 显著性水平
const DATA_SEED = 20261122; // 虚构数据生成种子（固定，不随分析种子变化）

/* 片区类型图：R 住宅 / C 商业 / I 工业 / . 非服务区（水域、绿地，不是片区） */
const TYPE_MAP = [
  'R.RRRRRR',
  '..RRRIIR',
  'RRCCRIIR',
  'RRCCRRRR',
  'RRCCRRRR',
  'RIIRRRRR',
  'RIIRRRR.',
  'RRRRRR.R',
];

const TYPE_NAME = { R: '住宅', C: '商业', I: '工业' };
const TYPE_BASE = { R: 80, C: 120, I: 95 }; // 各类型基线延迟(ms)：商业区天然更高

/* 指定若干片区为“无服务记录”与“覆盖不足”（虚构场景的一部分） */
const NO_DATA_IDS = new Set(['D7', 'F2']);
const LOW_COV = { B5: { n: 12, dayPool: 5 }, E7: { n: 18, dayPool: 6 }, G4: { n: 9, dayPool: 4 } };
/* 空间离群注入：B3 被高延迟邻居包围却自身很低（低高），G3 被低延迟邻居包围却自身很高（高低） */
const ANOMALY = { B3: -90, G3: 120 };

/* ---------------- 可重复随机数 ---------------- */
function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function makeGaussian(rng) {
  let spare = null;
  return function () {
    if (spare !== null) { const v = spare; spare = null; return v; }
    let u = 0, v = 0;
    while (u === 0) u = rng();
    v = rng();
    const r = Math.sqrt(-2 * Math.log(u)), th = 2 * Math.PI * v;
    spare = r * Math.sin(th);
    return r * Math.cos(th);
  };
}

/* ---------------- 虚构城市与数据生成 ---------------- */
function cellId(r, c) { return String.fromCharCode(65 + r) + (c + 1); }

function buildCity() {
  const districts = [];
  for (let r = 0; r < GRID; r++) {
    for (let c = 0; c < GRID; c++) {
      const t = TYPE_MAP[r][c];
      districts.push({
        idx: r * GRID + c, id: cellId(r, c), r, c,
        service: t !== '.', type: t === '.' ? null : t,
      });
    }
  }
  return districts;
}

/* 延迟场：两个高延迟热点 + 一个低延迟洼地，制造空间聚集 */
function hotspotField(r, c) {
  const g = (cr, cc, s) => Math.exp(-((r - cr) ** 2 + (c - cc) ** 2) / (2 * s * s));
  return 100 * g(1.5, 2.0, 1.6) + 65 * g(4.5, 3.5, 1.4) - 50 * g(6.0, 1.5, 1.3);
}

function genDataset() {
  const districts = buildCity();
  const rng = mulberry32(DATA_SEED);
  const gauss = makeGaussian(rng);
  let recId = 0;
  for (const d of districts) {
    d.records = [];
    if (!d.service) continue;
    if (NO_DATA_IDS.has(d.id)) continue; // 无服务记录
    const mean = Math.max(20, TYPE_BASE[d.type] + hotspotField(d.r, d.c) + (rng() - 0.5) * 24 + (ANOMALY[d.id] || 0));
    d.trueMean = mean;
    const low = LOW_COV[d.id];
    const n = low ? low.n : 30 + Math.floor(rng() * 55);
    for (let k = 0; k < n; k++) {
      const day = low ? Math.floor(rng() * low.dayPool) : Math.floor(rng() * DAYS);
      const latency = Math.max(5, mean * Math.exp(0.3 * gauss() - 0.045));
      d.records.push({
        rid: ++recId, day, hour: Math.floor(rng() * 24),
        latency: Math.round(latency * 10) / 10,
      });
    }
    d.records.sort((a, b) => a.day - b.day || a.hour - b.hour);
  }
  return districts;
}

/* 每片区：样本量、均值、有效覆盖 */
function districtStats(d) {
  const n = d.records.length;
  if (n === 0) return { n: 0, mean: null, cov: 0 };
  let s = 0, s2 = 0;
  const days = new Set();
  for (const r of d.records) { s += r.latency; days.add(r.day); }
  const mean = s / n;
  for (const r of d.records) s2 += (r.latency - mean) ** 2;
  return { n, mean, cov: days.size / DAYS, sd: Math.sqrt(s2 / Math.max(1, n - 1)) };
}

/* ---------------- 邻接 ---------------- */
function buildNeighbors(rule) {
  const offs = rule === 'queen'
    ? [[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]]
    : [[-1,0],[0,-1],[0,1],[1,0]];
  const map = new Map();
  for (const d of DISTRICTS) {
    if (!d.service) continue;
    const nb = [];
    for (const [dr, dc] of offs) {
      const r = d.r + dr, c = d.c + dc;
      if (r < 0 || r >= GRID || c < 0 || c >= GRID) continue;
      const nd = DISTRICTS[r * GRID + c];
      if (nd.service) nb.push(nd.idx);
    }
    map.set(d.idx, nb);
  }
  return map;
}

/* ---------------- 统计核心 ---------------- */
function meanOf(a) { let s = 0; for (const v of a) s += v; return s / a.length; }

/* 全局 Moran's I（行标准化权重） */
function moranGlobal(z, nbrsArr) {
  const m = z.length;
  let num = 0, sw = 0, szz = 0;
  for (let i = 0; i < m; i++) {
    szz += z[i] * z[i];
    const nb = nbrsArr[i], w = 1 / nb.length;
    for (const j of nb) { num += w * z[i] * z[j]; sw += w; }
  }
  return (m / sw) * (num / szz);
}

function fisherYatesPartial(arr, k, rng) {
  // 原地部分洗牌：洗完后前 k 个为不放回抽样结果
  const n = arr.length;
  for (let i = 0; i < Math.min(k, n - 1); i++) {
    const j = i + Math.floor(rng() * (n - i));
    const t = arr[i]; arr[i] = arr[j]; arr[j] = t;
  }
  return arr;
}

/* Benjamini–Hochberg 校正 */
function benjaminiHochberg(p) {
  const m = p.length;
  const order = p.map((_, i) => i).sort((a, b) => p[a] - p[b]);
  const q = new Array(m);
  let prev = 1;
  for (let k = m; k >= 1; k--) {
    const i = order[k - 1];
    prev = Math.min(prev, p[i] * m / k);
    q[i] = prev;
  }
  return q;
}

/* 主分析：返回全局 + 局部结果 */
function computeAnalysis(rule, stratified, nPerm, seed) {
  const nbrMap = buildNeighbors(rule);
  const reason = new Array(GRID * GRID).fill(null); // 排除原因
  const cand = [];
  for (const d of DISTRICTS) {
    if (!d.service) { reason[d.idx] = '非服务区'; continue; }
    const st = STATS[d.idx];
    if (st.n === 0) { reason[d.idx] = '无服务记录'; continue; }
    if (st.n < MIN_N || st.cov < MIN_COV) { reason[d.idx] = '覆盖不足'; continue; }
    cand.push(d.idx);
  }
  const candSet = new Set(cand);
  const nbrOf = {};
  for (const i of cand) nbrOf[i] = nbrMap.get(i).filter(j => candSet.has(j));
  const E = [];
  for (const i of cand) {
    if (nbrOf[i].length === 0) reason[i] = '无有效邻居';
    else E.push(i);
  }

  // 位置数组与 z 分数
  const pos = new Map(E.map((idx, k) => [idx, k]));
  const nbrsArr = E.map(i => nbrOf[i].map(j => pos.get(j)));
  const x = E.map(i => STATS[i].mean);
  const xbar = meanOf(x);
  const sd = Math.sqrt(meanOf(x.map(v => (v - xbar) ** 2)));
  const z = x.map(v => (v - xbar) / sd);
  const m = E.length;

  // 全局 Moran's I
  const Iobs = moranGlobal(z, nbrsArr);
  const rng = mulberry32(seed);
  const groups = stratified
    ? (() => { const g = {}; E.forEach((idx, k) => { (g[DISTRICTS[idx].type] ||= []).push(k); }); return Object.values(g); })()
    : [E.map((_, k) => k)];
  let cntG = 0;
  const zp = z.slice();
  for (let t = 0; t < nPerm; t++) {
    for (const g of groups) {
      for (let a = g.length - 1; a > 0; a--) {
        const b = Math.floor(rng() * (a + 1));
        const i1 = g[a], i2 = g[b];
        const tmp = zp[i1]; zp[i1] = zp[i2]; zp[i2] = tmp;
      }
    }
    if (moranGlobal(zp, nbrsArr) >= Iobs) cntG++;
  }
  const pGlobal = (1 + cntG) / (1 + nPerm);

  // 局部 LISA（条件随机化：固定 z_i，置换其余）
  const lisa = new Array(m);
  for (let k = 0; k < m; k++) {
    const nb = nbrsArr[k], kn = nb.length;
    let lag = 0;
    for (const j of nb) lag += z[j];
    lag /= kn;
    const Ii = z[k] * lag;
    // 置换池：同类型（可选）且排除自身
    const pool = [];
    for (let j = 0; j < m; j++) {
      if (j === k) continue;
      if (stratified && DISTRICTS[E[j]].type !== DISTRICTS[E[k]].type) continue;
      pool.push(z[j]);
    }
    let cnt = 0;
    const withReplacement = pool.length < kn; // 置换池小于邻居数时退化为有放回抽样
    for (let t = 0; t < nPerm; t++) {
      let lg = 0;
      if (withReplacement) {
        for (let j = 0; j < kn; j++) lg += pool[Math.floor(rng() * pool.length)];
      } else {
        const buf = pool.slice(); // 每次从原始池拷贝，保证不放回抽样的均匀性
        fisherYatesPartial(buf, kn, rng);
        for (let j = 0; j < kn; j++) lg += buf[j];
      }
      const It = z[k] * (lg / kn);
      if (Ii >= 0 ? It >= Ii : It <= Ii) cnt++;
    }
    lisa[k] = { lag, I: Ii, p: (1 + cnt) / (1 + nPerm), poolSize: pool.length };
  }
  const q = benjaminiHochberg(lisa.map(l => l.p));
  for (let k = 0; k < m; k++) {
    const L = lisa[k];
    L.q = q[k];
    L.sig = q[k] <= ALPHA;
    if (!L.sig || z[k] === 0 || L.lag === 0) L.cls = 'ns';
    else if (z[k] > 0 && L.lag > 0) L.cls = 'HH';
    else if (z[k] < 0 && L.lag < 0) L.cls = 'LL';
    else if (z[k] > 0 && L.lag < 0) L.cls = 'HL';
    else L.cls = 'LH';
  }

  return {
    rule, stratified, nPerm, seed, E, reason, nbrOf, nbrsArr,
    xbar, sd, z, Iobs, pGlobal, cntG, lisa,
    expectedI: -1 / (m - 1),
    minP: 1 / (nPerm + 1),
    rawSig: lisa.filter(l => l.p <= ALPHA).length,
    bhSig: lisa.filter(l => l.sig).length,
    counts: {
      HH: lisa.filter(l => l.cls === 'HH').length,
      LL: lisa.filter(l => l.cls === 'LL').length,
      HL: lisa.filter(l => l.cls === 'HL').length,
      LH: lisa.filter(l => l.cls === 'LH').length,
    },
  };
}

/* ---------------- 成对样本：同一多重集，两种空间排列 ---------------- */
function pairedSample(rule, nPerm, seed) {
  const nbrMap = buildNeighbors(rule);
  const P = DISTRICTS.filter(d => d.service && nbrMap.get(d.idx).length > 0).map(d => d.idx);
  const pos = new Map(P.map((idx, k) => [idx, k]));
  const nbrsArr = P.map(i => nbrMap.get(i).map(j => pos.get(j)));
  // 确定性取值场（聚集）
  const val = i => {
    const d = DISTRICTS[i];
    const g = (cr, cc, s) => Math.exp(-((d.r - cr) ** 2 + (d.c - cc) ** 2) / (2 * s * s));
    return Math.round((60 + 90 * g(2, 2, 1.8) + 55 * g(5.5, 5, 1.5) + 3 * Math.sin(d.r * 2.1 + d.c * 1.3)) * 100) / 100;
  };
  const vA = P.map(val);
  // B：同样的值，黑白格交错（高值给黑格、低值给白格）
  const blacks = P.filter(i => (DISTRICTS[i].r + DISTRICTS[i].c) % 2 === 0);
  const whites = P.filter(i => (DISTRICTS[i].r + DISTRICTS[i].c) % 2 === 1);
  const sorted = vA.slice().sort((a, b) => b - a);
  const vB = new Array(P.length);
  blacks.forEach((idx, k) => { vB[pos.get(idx)] = sorted[k]; });
  whites.forEach((idx, k) => { vB[pos.get(idx)] = sorted[blacks.length + k]; });

  const rng = mulberry32(seed);
  const run = v => {
    const xbar = meanOf(v);
    const sd = Math.sqrt(meanOf(v.map(t => (t - xbar) ** 2)));
    const z = v.map(t => (t - xbar) / sd);
    const I = moranGlobal(z, nbrsArr);
    let cnt = 0;
    const zp = z.slice();
    for (let t = 0; t < nPerm; t++) {
      for (let a = zp.length - 1; a > 0; a--) {
        const b = Math.floor(rng() * (a + 1));
        const tmp = zp[a]; zp[a] = zp[b]; zp[b] = tmp;
      }
      if (moranGlobal(zp, nbrsArr) >= I) cnt++;
    }
    return { I, p: (1 + cnt) / (1 + nPerm), mean: xbar, sd, z };
  };
  return { P, nbrsArr, vA, vB, A: run(vA), B: run(vB) };
}

/* ================= 渲染 ================= */
const DISTRICTS = genDataset();
const STATS = DISTRICTS.map(districtStats);

const state = { rule: 'rook', stratified: false, nPerm: 4999, seed: 20261122, selected: null };
let ANAL = null;   // 主分析结果
let PAIR = null;   // 成对样本结果

const $ = id => document.getElementById(id);
function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}
const fmt = (v, d = 1) => v === null || v === undefined ? '—' : Number(v).toFixed(d);
const fmtP = p => p <= ANAL.minP ? `≤ ${ANAL.minP.toFixed(4)}（下限）` : p.toFixed(4);

/* 顺序色阶（8 档，CSS 变量驱动，深色模式自动翻转锚点） */
const SEQ_VARS = ['--seq-1','--seq-2','--seq-3','--seq-4','--seq-5','--seq-6','--seq-7','--seq-8'];
function seqColor(t) {
  const i = Math.max(0, Math.min(7, Math.floor(t * 7.999)));
  return `var(${SEQ_VARS[i]})`;
}

/* ---------- 提示框 ---------- */
const tooltip = () => $('tooltip');
function showTip(html, x, y) {
  const tp = tooltip();
  tp.innerHTML = '';
  for (const [txt, cls] of html) tp.appendChild(el('div', cls, txt));
  tp.hidden = false;
  const pad = 12;
  const r = tp.getBoundingClientRect();
  tp.style.left = Math.min(x + pad, innerWidth - r.width - 8) + 'px';
  tp.style.top = Math.min(y + pad, innerHeight - r.height - 8) + 'px';
}
function hideTip() { tooltip().hidden = true; }

/* ---------- 地图构建 ---------- */
function makeCell(d) {
  const b = el('button', 'cell');
  b.dataset.idx = d.idx;
  b.setAttribute('aria-label', `片区 ${d.id}`);
  if (d.service) {
    b.appendChild(el('span', 'badge', TYPE_NAME[d.type][0]));
  } else {
    b.classList.add('nonservice');
    b.disabled = true;
  }
  return b;
}

function renderMeanMap() {
  const map = $('meanMap');
  map.innerHTML = '';
  const vals = ANAL.E.map(i => STATS[i].mean);
  const lo = Math.min(...vals), hi = Math.max(...vals);
  const inE = new Set(ANAL.E);
  for (const d of DISTRICTS) {
    const c = makeCell(d);
    if (d.service) {
      if (inE.has(d.idx)) {
        const t = (STATS[d.idx].mean - lo) / (hi - lo || 1);
        c.style.background = seqColor(t);
      } else {
        c.classList.add('excluded');
      }
      const tipLines = () => {
        const st = STATS[d.idx];
        return [
          [`片区 ${d.id} · ${TYPE_NAME[d.type]}`, 'tv'],
          [st.n ? `均值 ${fmt(st.mean)} ms · n=${st.n} · 覆盖 ${(st.cov * 100).toFixed(0)}%` : '无服务记录', 'tl'],
          [ANAL.reason[d.idx] ? `不参与计算：${ANAL.reason[d.idx]}` : '参与计算', 'tl'],
        ];
      };
      c.addEventListener('pointermove', ev => showTip(tipLines(), ev.clientX, ev.clientY));
      c.addEventListener('pointerleave', hideTip);
      c.addEventListener('focus', () => {
        const r = c.getBoundingClientRect();
        showTip(tipLines(), r.left + r.width / 2, r.bottom);
      });
      c.addEventListener('blur', hideTip);
      c.addEventListener('click', () => selectCell(d.idx));
    }
    if (state.selected === d.idx) c.classList.add('selected');
    map.appendChild(c);
  }
  // 图例 ramp
  const ramp = $('meanRamp');
  ramp.innerHTML = '';
  for (let i = 0; i < 8; i++) {
    const s = el('span');
    s.style.background = `var(${SEQ_VARS[i]})`;
    ramp.appendChild(s);
  }
  const lab = $('meanRampLabels');
  lab.innerHTML = '';
  lab.appendChild(el('span', null, `${fmt(lo)} ms`));
  lab.appendChild(el('span', null, `${fmt(hi)} ms`));
}

const CLS_LABEL = { HH: '高高聚集', LL: '低低聚集', HL: '高低离群', LH: '低高离群', ns: '不显著' };
const CLS_VAR = { HH: '--hh', LL: '--ll', HL: '--hl', LH: '--lh' };

function renderLisaMap() {
  const map = $('lisaMap');
  map.innerHTML = '';
  const pos = new Map(ANAL.E.map((idx, k) => [idx, k]));
  for (const d of DISTRICTS) {
    const c = makeCell(d);
    if (d.service) {
      const k = pos.get(d.idx);
      if (k === undefined) {
        c.classList.add('excluded');
      } else {
        const L = ANAL.lisa[k];
        if (L.cls === 'ns') {
          c.style.background = 'var(--ns)';
        } else {
          c.style.background = `var(${CLS_VAR[L.cls]})`;
          const lab = el('span', 'cls', L.cls);
          lab.style.color = L.cls === 'HL' ? '#0b0b0b' : '#fff';
          c.appendChild(lab);
        }
      }
      const tipLines = () => {
        const k2 = pos.get(d.idx);
        const lines = [[`片区 ${d.id} · ${TYPE_NAME[d.type]}`, 'tv']];
        if (k2 === undefined) lines.push([`不参与计算：${ANAL.reason[d.idx]}`, 'tl']);
        else {
          const L = ANAL.lisa[k2];
          lines.push([`${CLS_LABEL[L.cls]} · Iᵢ=${fmt(L.I, 3)} · q=${L.q.toFixed(4)}`, 'tl']);
        }
        return lines;
      };
      c.addEventListener('pointermove', ev => showTip(tipLines(), ev.clientX, ev.clientY));
      c.addEventListener('pointerleave', hideTip);
      c.addEventListener('focus', () => {
        const r = c.getBoundingClientRect();
        showTip(tipLines(), r.left + r.width / 2, r.bottom);
      });
      c.addEventListener('blur', hideTip);
      c.addEventListener('click', () => selectCell(d.idx));
    }
    if (state.selected === d.idx) c.classList.add('selected');
    map.appendChild(c);
  }
  // 图例
  const lg = $('lisaLegend');
  lg.innerHTML = '';
  const items = [
    ['HH 高高聚集', '--hh'], ['LL 低低聚集', '--ll'],
    ['HL 高低离群', '--hl'], ['LH 低高离群', '--lh'],
  ];
  for (const [txt, v] of items) {
    const s = el('span');
    const sw = el('span', 'sw'); sw.style.background = `var(${v})`;
    s.appendChild(sw); s.appendChild(document.createTextNode(txt));
    lg.appendChild(s);
  }
  const ns = el('span');
  const sw2 = el('span', 'sw'); sw2.style.background = 'var(--ns)';
  ns.appendChild(sw2); ns.appendChild(document.createTextNode('不显著'));
  lg.appendChild(ns);
  const ex = el('span');
  ex.appendChild(el('span', 'sw hatch'));
  ex.appendChild(document.createTextNode('不参与计算'));
  lg.appendChild(ex);
  $('lisaNote').textContent =
    `原始 p ≤ 0.05 的片区有 ${ANAL.rawSig} 个，经 BH 多重比较校正后显著的为 ${ANAL.bhSig} 个 —— 未校正的原始小概率不直接标红。`;
}

/* ---------- 统计磁贴 ---------- */
function renderTiles() {
  const t = $('tiles');
  t.innerHTML = '';
  const pTxt = ANAL.pGlobal <= ANAL.minP
    ? `≤ ${ANAL.minP.toFixed(4)}` : ANAL.pGlobal.toFixed(4);
  const tiles = [
    [fmt(ANAL.Iobs, 3), '全局 Moran’s I', `随机期望 E[I] = ${fmt(ANAL.expectedI, 3)}`],
    [pTxt, '置换 p 值（单侧·正向聚集）', `最小可分辨 p = 1/(${ANAL.nPerm}+1) = ${ANAL.minP.toFixed(4)}`],
    [`${ANAL.E.length} / ${GRID * GRID}`, '参与计算的片区', '其余片区单列，不当作零值'],
    [String(ANAL.counts.HH), '高高聚集片区（BH 校正后）', `另有低低聚集 ${ANAL.counts.LL} 个`],
    [String(ANAL.counts.HL + ANAL.counts.LH), '空间离群片区（高低+低高）', `HL ${ANAL.counts.HL} · LH ${ANAL.counts.LH}`],
    [`${ANAL.rawSig} → ${ANAL.bhSig}`, '原始 p≤0.05 → BH 校正后', '多重比较校正，避免满图标红'],
  ];
  for (const [v, l, s] of tiles) {
    const d = el('div', 'tile');
    d.appendChild(el('div', 'v', v));
    d.appendChild(el('div', 'l', l));
    d.appendChild(el('div', 's', s));
    t.appendChild(d);
  }
}

/* ---------- 片区详情 ---------- */
function selectCell(idx) {
  state.selected = idx;
  renderMeanMap();
  renderLisaMap();
  renderDetail();
}

function renderDetail() {
  const body = $('detailBody');
  const title = $('detailTitle');
  body.innerHTML = '';
  const idx = state.selected;
  if (idx === null) {
    title.textContent = '片区详情';
    body.appendChild(el('p', 'note', '点击上方任一地图格子，查看该片区的实际邻居、权重、参与计算的服务记录，以及它属于高高聚集 / 低低聚集 / 高低离群 / 低高离群中的哪一类。更换邻接规则后，此处内容同步变化。'));
    return;
  }
  const d = DISTRICTS[idx];
  const st = STATS[idx];
  title.textContent = `片区 ${d.id} · ${d.service ? TYPE_NAME[d.type] + '区' : '非服务区'}`;
  if (!d.service) {
    body.appendChild(el('p', 'note', '该格子为非服务区（水域 / 绿地），不是片区，不参与任何统计。'));
    return;
  }
  const pos = new Map(ANAL.E.map((i, k) => [i, k]));
  const k = pos.get(idx);

  const grid = el('div', 'detail-grid');
  const left = el('div'), right = el('div');
  grid.appendChild(left); grid.appendChild(right);
  body.appendChild(grid);

  // 左：基本量 + LISA 结果
  left.appendChild(el('h3', null, '片区统计量'));
  const kv = el('dl', 'kv');
  const addKV = (a, b) => { kv.appendChild(el('dt', null, a)); kv.appendChild(el('dd', null, b)); };
  addKV('样本量 n', String(st.n));
  addKV('平均延迟', st.n ? `${fmt(st.mean)} ms` : '—');
  addKV('有效覆盖', st.n ? `${(st.cov * 100).toFixed(0)}%（${Math.round(st.cov * DAYS)}/${DAYS} 天有记录）` : '—');
  if (k !== undefined) {
    addKV('z 分数', fmt(ANAL.z[k], 3));
    addKV('空间滞后（邻居加权均值）', fmt(ANAL.lisa[k].lag, 3));
    addKV('局部 Moran Iᵢ', fmt(ANAL.lisa[k].I, 3));
    addKV('原始 p', fmtP(ANAL.lisa[k].p));
    addKV('BH 校正 q', ANAL.lisa[k].q.toFixed(4));
  }
  left.appendChild(kv);

  if (k !== undefined) {
    const L = ANAL.lisa[k];
    const pill = el('span', `pill ${L.cls === 'ns' ? 'ns' : L.cls.toLowerCase()}`, CLS_LABEL[L.cls]);
    const p1 = el('p', 'note');
    p1.appendChild(document.createTextNode('分类（BH 校正后）：'));
    p1.appendChild(pill);
    left.appendChild(p1);
    const expl = {
      HH: '本片区延迟高于全局均值，邻居加权均值也偏高 → 高值被高值包围。',
      LL: '本片区延迟低于全局均值，邻居加权均值也偏低 → 低值被低值包围。',
      HL: '本片区延迟高于全局均值，但邻居加权均值偏低 → 高值被低值包围的空间离群。',
      LH: '本片区延迟低于全局均值，但邻居加权均值偏高 → 低值被高值包围的空间离群。',
      ns: '校正后不显著：在该邻接规则与置换检验下，不能断言它与邻居存在空间关联。',
    }[L.cls];
    left.appendChild(el('p', 'note', expl));
    if (state.stratified && L.poolSize < 10) {
      const w = el('p', 'note warn-text',
        `注意：同类型置换池仅 ${L.poolSize} 个片区，置换精度受限，解释 p 值需谨慎。`);
      left.appendChild(w);
    }
  } else {
    const w = el('p', 'note warn-text', `该片区不参与计算：${ANAL.reason[idx]}。它不是零值，从全局与局部统计中剔除并单列。`);
    left.appendChild(w);
  }

  // 左下：参与计算的服务记录
  left.appendChild(el('h3', null, `参与计算的服务记录（共 ${st.n} 条）`));
  if (st.n === 0) {
    left.appendChild(el('p', 'note', '该片区没有任何服务记录。'));
  } else {
    const tbl = el('table');
    const thead = el('thead');
    const hr = el('tr');
    for (const h of ['记录号', '第几天', '时段', '延迟 (ms)']) hr.appendChild(el('th', null, h));
    thead.appendChild(hr); tbl.appendChild(thead);
    const tb = el('tbody');
    for (const r of d.records.slice(0, 12)) {
      const tr = el('tr');
      tr.appendChild(el('td', null, `#${r.rid}`));
      tr.appendChild(el('td', null, `D${r.day + 1}`));
      tr.appendChild(el('td', null, `${String(r.hour).padStart(2, '0')}:00`));
      tr.appendChild(el('td', null, fmt(r.latency)));
      tb.appendChild(tr);
    }
    tbl.appendChild(tb);
    left.appendChild(tbl);
    if (st.n > 12) left.appendChild(el('p', 'note', `仅列出前 12 条，全部 ${st.n} 条均参与均值计算。`));
  }

  // 右：邻居与权重（随邻接规则变化）
  const ruleName = state.rule === 'rook' ? '共享边' : '共享角';
  right.appendChild(el('h3', null, `实际邻居与权重（${ruleName}规则）`));
  const svcNb = buildNeighbors(state.rule).get(idx) || [];
  if (svcNb.length === 0) {
    right.appendChild(el('p', 'note', `在${ruleName}规则下，该片区没有任何服务区邻居 → 列入“无有效邻居”，不参与空间统计。`));
  } else {
    const eSet = new Set(ANAL.E);
    const effNb = svcNb.filter(j => eSet.has(j));
    const w = 1 / Math.max(1, effNb.length);
    const tbl = el('table');
    const thead = el('thead');
    const hr = el('tr');
    for (const h of ['邻居', '类型', '权重 w', '均值 (ms)', 'z']) hr.appendChild(el('th', null, h));
    thead.appendChild(hr); tbl.appendChild(thead);
    const tb = el('tbody');
    for (const j of svcNb) {
      const nd = DISTRICTS[j];
      const tr = el('tr');
      tr.appendChild(el('td', null, nd.id));
      tr.appendChild(el('td', null, TYPE_NAME[nd.type]));
      if (eSet.has(j)) {
        const kj = pos.get(j);
        tr.appendChild(el('td', null, `1/${effNb.length} ≈ ${w.toFixed(3)}`));
        tr.appendChild(el('td', null, fmt(STATS[j].mean)));
        tr.appendChild(el('td', null, fmt(ANAL.z[kj], 2)));
      } else {
        tr.appendChild(el('td', null, '—'));
        tr.appendChild(el('td', null, STATS[j].n ? fmt(STATS[j].mean) : '—'));
        const td = el('td'); td.appendChild(el('span', 'pill ex', ANAL.reason[j] || '不参与'));
        tr.appendChild(td);
      }
      tb.appendChild(tr);
    }
    tbl.appendChild(tb);
    right.appendChild(tbl);
    right.appendChild(el('p', 'note',
      `行标准化权重：每个有效邻居权重 = 1 / 有效邻居数（此处为 ${effNb.length}）。被排除的邻居不分配权重、不当作零值。`));
  }
}

/* ---------- 排除清单 ---------- */
function renderExcluded() {
  const row = $('excludedRow');
  row.innerHTML = '';
  const groups = { '无服务记录': [], '覆盖不足': [], '无有效邻居': [] };
  for (const d of DISTRICTS) {
    const r = ANAL.reason[d.idx];
    if (r && groups[r]) groups[r].push(d);
  }
  const defs = [
    ['无服务记录', g => `共 ${g.length} 个：${g.map(d => d.id).join('、') || '无'}`,
     '这些片区一条服务记录都没有。缺失不是“延迟为 0”，它们从均值、z 分数与所有空间统计中剔除并单列。'],
    ['覆盖不足', g => `共 ${g.length} 个：${g.map(d => `${d.id}（n=${STATS[d.idx].n}，覆盖 ${(STATS[d.idx].cov * 100).toFixed(0)}%）`).join('、') || '无'}`,
     `参与计算要求 n ≥ ${MIN_N} 且有效覆盖 ≥ ${MIN_COV * 100}%（${DAYS} 天观察窗内有记录的天数占比）。样本太少时均值抽样误差大（SE ≈ s/√n），覆盖太低时样本代表不了整个观察窗；把这样的均值当作真值去做置换检验，会污染自身与邻居的统计量，因此不能下结论。`],
    ['无有效邻居', g => `共 ${g.length} 个：${g.map(d => d.id).join('、') || '无'}`,
     `在当前邻接规则（${state.rule === 'rook' ? '共享边' : '共享角'}）下，这些片区四周没有可参与计算的邻居，空间权重无法定义。它们不是零值：全局统计剔除，局部统计无定义。切换邻接规则后此名单会变化。`],
  ];
  for (const [name, countFn, expl] of defs) {
    const box = el('div', 'excl-box');
    box.appendChild(el('h3', null, name));
    box.appendChild(el('div', 'ids', countFn(groups[name])));
    box.appendChild(el('p', null, expl));
    row.appendChild(box);
  }
}

/* ---------- 成对样本 ---------- */
function renderPaired() {
  const renderMini = (mapId, P, v) => {
    const map = $(mapId);
    map.innerHTML = '';
    const lo = Math.min(...v), hi = Math.max(...v);
    const valOf = new Map(P.map((idx, k) => [idx, v[k]]));
    for (const d of DISTRICTS) {
      const c = el('div', 'cell');
      if (!d.service) c.classList.add('nonservice');
      else if (!valOf.has(d.idx)) c.classList.add('excluded');
      else {
        c.style.background = seqColor((valOf.get(d.idx) - lo) / (hi - lo || 1));
        c.title = `${d.id}: ${fmt(valOf.get(d.idx), 1)} ms`;
      }
      c.style.cursor = 'default';
      map.appendChild(c);
    }
  };
  renderMini('pairMapA', PAIR.P, PAIR.vA);
  renderMini('pairMapB', PAIR.P, PAIR.vB);

  const box = $('pairedStats');
  box.innerHTML = '';
  const tbl = el('table');
  const thead = el('thead');
  const hr = el('tr');
  for (const h of ['', 'A · 聚集排列', 'B · 交错排列']) hr.appendChild(el('th', null, h));
  thead.appendChild(hr); tbl.appendChild(thead);
  const tb = el('tbody');
  const rows = [
    ['全局均值 (ms)', fmt(PAIR.A.mean, 2), fmt(PAIR.B.mean, 2)],
    ['标准差 (ms)', fmt(PAIR.A.sd, 2), fmt(PAIR.B.sd, 2)],
    ['Moran’s I', fmt(PAIR.A.I, 3), fmt(PAIR.B.I, 3)],
    ['置换 p（单侧）',
      PAIR.A.p <= PAIR.minP ? `≤ ${PAIR.minP.toFixed(4)}` : PAIR.A.p.toFixed(4),
      PAIR.B.p <= PAIR.minP ? `≤ ${PAIR.minP.toFixed(4)}` : PAIR.B.p.toFixed(4)],
  ];
  for (const [a, b, c] of rows) {
    const tr = el('tr');
    tr.appendChild(el('td', null, a));
    tr.appendChild(el('td', null, b));
    tr.appendChild(el('td', null, c));
    tb.appendChild(tr);
  }
  tbl.appendChild(tb);
  box.appendChild(tbl);
  const concl = el('p', 'note');
  concl.appendChild(document.createTextNode(
    `两组均值完全相同（${fmt(PAIR.A.mean, 2)} ms），但 Moran’s I 分别为 ${fmt(PAIR.A.I, 3)} 与 ${fmt(PAIR.B.I, 3)}：` +
    `A 呈显著空间聚集，B 不聚集（甚至趋于交错离散）。空间结论由“值放在哪里”决定，不能仅由平均值决定。`));
  box.appendChild(concl);
  box.appendChild(el('p', 'note', `成对样本使用与主分析相同的邻接规则（${state.rule === 'rook' ? '共享边' : '共享角'}）、置换次数与随机种子；无邻居格子同样剔除。`));
}

/* ---------- 方法说明 ---------- */
function renderMethod() {
  const dl = $('methodList');
  dl.innerHTML = '';
  const items = [
    ['随机种子与可重复性', `当前分析种子 = ${state.seed}。同一种子 + 同一置换次数 + 同一邻接规则会得到完全相同的置换序列与 p 值；改种子只改变随机序列，不改变数据。虚构数据集本身由固定生成种子 ${DATA_SEED} 产生。`],
    ['置换次数与最小精度', `当前置换 ${state.nPerm} 次，可估计的最小 p = 1/(${state.nPerm}+1) ≈ ${ANAL.minP.toFixed(4)}。当观察统计量落在所有置换之外时，只能报告 p ≤ ${ANAL.minP.toFixed(4)}，更小的真实 p 需要更多置换才能分辨。`],
    ['多重比较校正', `局部检验同时对 ${ANAL.E.length} 个片区做假设检验。即使全无空间关联，也期望约 ${(ANAL.E.length * ALPHA).toFixed(1)} 个片区仅凭运气得到 p ≤ 0.05。因此用 Benjamini–Hochberg 控制错误发现率：本次原始 p ≤ 0.05 的有 ${ANAL.rawSig} 个，校正后显著 ${ANAL.bhSig} 个，只有校正后的才在地图上着色。`],
    ['同类型片区内置换', `商业区（基线 120ms）与住宅区（80ms）、工业区（95ms）的延迟分布本来就不同。自由置换会把“类型差异”混进“空间关联”；同类型内置换只在同类型片区之间打乱取值，回答的是“在同类型中，相邻的是否更相似”。`],
    ['为什么覆盖不足不能下结论', `参与计算要求 n ≥ ${MIN_N} 且有效覆盖 ≥ ${MIN_COV * 100}%。n 太小时，均值的抽样标准误 SE ≈ s/√n 很大（例如 n=9 时 SE 是 n=81 时的 3 倍）；覆盖太低时，记录只来自 ${DAYS} 天中的少数几天，无法代表整个观察窗。把这样的均值当作真值进入置换检验，既骗自己也污染邻居——所以覆盖不足的片区一律剔除并单列，而不是当作零值或正常值。`],
    ['无邻居与无记录片区', '无服务记录的片区没有任何可统计的观测；无有效邻居的片区无法定义空间权重。两者都不是“延迟为 0”，全部从统计中剔除并单列，且名单随邻接规则变化。'],
    ['统计量定义', `全局 Moran’s I 使用行标准化权重，随机期望 E[I] = -1/(m-1) = ${fmt(ANAL.expectedI, 4)}（m = 参与片区数）。局部 Iᵢ = zᵢ × 邻居 z 的加权均值；显著性用条件随机化（固定 zᵢ，置换其余）估计，p 为 (1+超过观察值的置换数)/(1+置换次数)。`],
  ];
  for (const [t, d] of items) {
    dl.appendChild(el('dt', null, t));
    dl.appendChild(el('dd', null, d));
  }
}

/* ---------- 全量表格 ---------- */
function renderTable() {
  const tbl = $('fullTable');
  tbl.innerHTML = '';
  const thead = el('thead');
  const hr = el('tr');
  for (const h of ['片区', '类型', 'n', '均值 (ms)', '有效覆盖', '状态', 'z', 'Iᵢ', '原始 p', 'BH q', '分类'])
    hr.appendChild(el('th', null, h));
  thead.appendChild(hr); tbl.appendChild(thead);
  const tb = el('tbody');
  const pos = new Map(ANAL.E.map((i, k) => [i, k]));
  for (const d of DISTRICTS) {
    const st = STATS[d.idx];
    const tr = el('tr');
    tr.appendChild(el('td', null, d.id));
    tr.appendChild(el('td', null, d.service ? TYPE_NAME[d.type] : '—'));
    tr.appendChild(el('td', null, d.service ? String(st.n) : '—'));
    tr.appendChild(el('td', null, d.service && st.n ? fmt(st.mean) : '—'));
    tr.appendChild(el('td', null, d.service && st.n ? `${(st.cov * 100).toFixed(0)}%` : '—'));
    const status = !d.service ? '非服务区' : (ANAL.reason[d.idx] || '参与计算');
    const tdS = el('td');
    tdS.appendChild(el('span', `pill ${status === '参与计算' ? 'ns' : 'ex'}`, status));
    tr.appendChild(tdS);
    const k = pos.get(d.idx);
    if (k !== undefined) {
      const L = ANAL.lisa[k];
      tr.appendChild(el('td', null, fmt(ANAL.z[k], 2)));
      tr.appendChild(el('td', null, fmt(L.I, 3)));
      tr.appendChild(el('td', null, L.p <= ANAL.minP ? `≤${ANAL.minP.toFixed(4)}` : L.p.toFixed(4)));
      tr.appendChild(el('td', null, L.q.toFixed(4)));
      const tdC = el('td');
      tdC.appendChild(el('span', `pill ${L.cls === 'ns' ? 'ns' : L.cls.toLowerCase()}`, CLS_LABEL[L.cls]));
      tr.appendChild(tdC);
    } else {
      for (let i = 0; i < 5; i++) tr.appendChild(el('td', null, '—'));
    }
    tb.appendChild(tr);
  }
  tbl.appendChild(tb);
}

/* ---------------- 主流程 ---------------- */
function recompute() {
  ANAL = computeAnalysis(state.rule, state.stratified, state.nPerm, state.seed);
  PAIR = pairedSample(state.rule, state.nPerm, state.seed);
  PAIR.minP = ANAL.minP;
  renderTiles();
  renderMeanMap();
  renderLisaMap();
  renderDetail();
  renderExcluded();
  renderPaired();
  renderMethod();
  renderTable();
}

function init() {
  document.querySelectorAll('input[name="adj"]').forEach(r =>
    r.addEventListener('change', () => { state.rule = r.value; recompute(); }));
  document.querySelectorAll('input[name="perm"]').forEach(r =>
    r.addEventListener('change', () => { state.stratified = r.value === 'strat'; recompute(); }));
  $('nPerm').addEventListener('change', e => { state.nPerm = parseInt(e.target.value, 10); recompute(); });
  $('seed').addEventListener('change', e => {
    const v = parseInt(e.target.value, 10);
    state.seed = Number.isFinite(v) ? v : 0;
    recompute();
  });
  $('themeToggle').addEventListener('click', () => {
    const root = document.documentElement;
    const dark = root.getAttribute('data-theme') === 'dark';
    root.setAttribute('data-theme', dark ? 'light' : 'dark');
    $('themeToggle').textContent = dark ? '切换深色' : '切换浅色';
  });
  recompute();
}

if (typeof document !== 'undefined') init();

/* 供测试引用的纯函数出口 */
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    mulberry32, buildCity, genDataset, districtStats, buildNeighbors,
    moranGlobal, benjaminiHochberg, computeAnalysis, pairedSample,
    DISTRICTS, STATS, MIN_N, MIN_COV, DAYS, GRID,
  };
}
