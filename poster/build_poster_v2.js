// A0 portrait poster, two-column version
const pptxgen = require("pptxgenjs");
const path = require("path");
const { execSync } = require("child_process");
const { applyTheme } = require("/root/.claude/skills/synced/f6ad0e85-3217-498c-aa06-38e6c8806359_0ff8945c-f50f-42d2-854f-a9da0a03ddaa/pptx/scripts/apply_theme.js");

const IMG = (f) => path.join(__dirname, "img", f);
const SIZE = JSON.parse(execSync(`python3 -c "
import json,glob,os
from PIL import Image
print(json.dumps({os.path.basename(f): Image.open(f).size for f in glob.glob('${path.join(__dirname, "img")}/*.png')}))"`).toString());

const THEME = {
  name: "Rn222 Poster",
  headFontFace: "Calibri", bodyFontFace: "Calibri",
  colors: {
    dk1: "1B2A4A", lt1: "FFFFFF", dk2: "1B2A4A", lt2: "F2F5F8",
    accent1: "1F8A8A", accent2: "D9631E", accent3: "6B7385", accent4: "C3CFDF",
    accent5: "E8B23A", accent6: "D5ECEC", hlink: "1F8A8A", folHlink: "6B7385",
  },
};

const pres = new pptxgen();
pres.defineLayout({ name: "A0P", width: 33.11, height: 46.81 });
pres.layout = "A0P";
pres.theme = { headFontFace: "Calibri", bodyFontFace: "Calibri" };
pres.title = "Measuring the Radon-222 Background in the XENONnT Dark Matter Detector";
pres.author = "Jie Tang";
const C = pres.SchemeColor;
pres.addSection({ title: "Poster" });
const s = pres.addSlide({ sectionTitle: "Poster" });
s.background = { color: "FFFFFF" };

// ---------------------------------------------------------------- geometry & type
const PW = 33.11, PH = 46.81, M = 0.9, GAP = 0.8;
const HDR_H = 5.9, FTR_H = 1.9;
const COLW = (PW - 2 * M - GAP) / 2;
const colX = [M, M + COLW + GAP];
const PAD = 0.5, INW = COLW - 2 * PAD;
const BODY = 30, CAP = 22, SUB = 34, HEAD = 52;
const PT = 1 / 72;
const ROWGAP = 0.45;

function nLines(text, size, width) {
  const cpl = Math.max(8, Math.floor(width / (size * PT * 0.43)));
  return text.split("\n").reduce((n, p) => n + Math.max(1, Math.ceil(p.length / cpl)), 0);
}
const lh = (size, sp = 1.17) => size * PT * sp;
const t = (text, o = {}) => ({ text, options: o });
const sup = (x) => t(x, { superscript: true });
const plain = (runs) => (typeof runs === "string" ? runs : runs.map((r) => r.text).join(""));

// ---------------------------------------------------------------- blocks: (w) => {h, draw(x,y)}
const para = (runs, o = {}) => (w) => {
  const size = o.size || BODY;
  const h = nLines(plain(runs), size, w) * lh(size) + 0.1;
  return { h: h + (o.after ?? 0.3), draw: (x, y) => s.addText(runs, { x, y, w, h, fontSize: size, color: o.color || C.text1,
    italic: o.italic, margin: 0, valign: "top", isTextBox: true, lineSpacingMultiple: 1.05 }) };
};
const bullets = (items, o = {}) => (w) => {
  const size = o.size || BODY;
  const runs = [];
  let nl = 0;
  items.forEach((it, i) => {
    const parts = typeof it === "string" ? [t(it)] : it;
    parts.forEach((p, j) => {
      const op = Object.assign({}, p.options);
      if (j === 0) op.bullet = { indent: 30 };
      if (j === parts.length - 1 && i < items.length - 1) op.breakLine = true;
      op.paraSpaceAfter = 10;
      runs.push(t(p.text, op));
    });
    nl += nLines(plain(parts), size, w - 0.45);
  });
  const h = nl * lh(size) + items.length * 10 * PT + 0.1;
  return { h: h + (o.after ?? 0.3), draw: (x, y) => s.addText(runs, { x, y, w, h, fontSize: size, color: C.text1,
    margin: 0, valign: "top", isTextBox: true, lineSpacingMultiple: 1.05 }) };
};
const sub = (text, color) => (w) => {
  const h = lh(SUB) + 0.05;
  return { h: h + 0.15, draw: (x, y) => s.addText(text, { x, y, w, h, fontSize: SUB, bold: true, color: color || C.accent1,
    margin: 0, valign: "top", isTextBox: true }) };
};
const image = (file, o = {}) => (w) => {
  const [pw, ph] = SIZE[file];
  const iw = o.frac ? w * o.frac : w;
  const h = iw * ph / pw;
  return { h: h + (o.after ?? 0.12), draw: (x, y) => s.addImage({ path: IMG(file), x: x + (w - iw) / 2, y, w: iw, h,
    altText: o.alt || file, objectName: o.name || file }) };
};
const caption = (text) => (w) => {
  const h = nLines(text, CAP, w) * lh(CAP, 1.12) + 0.06;
  return { h: h + 0.3, draw: (x, y) => s.addText(text, { x, y, w, h, fontSize: CAP, italic: true, color: C.accent3,
    margin: 0, valign: "top", isTextBox: true }) };
};
const placeholder = (label, hh) => (w) => ({ h: hh + 0.12, draw: (x, y) => {
  s.addShape(pres.shapes.RECTANGLE, { x, y, w, h: hh, fill: { color: C.background1 }, line: { color: C.accent3, width: 2, dashType: "dash" }, objectName: label });
  s.addText(label, { x, y, w, h: hh, fontSize: 26, color: C.accent3, align: "center", valign: "middle", isTextBox: true });
} });
// side-by-side: cols = [[widthFraction, [blocks...]], ...]
const row = (cols, o = {}) => (w) => {
  const gap = o.gap ?? 0.45;
  const widths = cols.map(([f]) => (w - gap * (cols.length - 1)) * f);
  const stacks = cols.map(([, blocks], i) => blocks.map((b) => b(widths[i])));
  const hs = stacks.map((st) => st.reduce((a, b) => a + b.h, 0));
  const H = Math.max(...hs);
  return { h: H + (o.after ?? 0.15), draw: (x, y) => {
    let xx = x;
    stacks.forEach((st, i) => {
      let yy = y + (o.vcenter ? (H - hs[i]) / 2 : 0);
      st.forEach((b) => { b.draw(xx, yy); yy += b.h; });
      xx += widths[i] + gap;
    });
  } };
};

// ---------------------------------------------------------------- section frame
const cursor = [HDR_H + 0.6, HDR_H + 0.6];
function section(col, num, title, blockFns) {
  const headH = 1.0;
  const blocks = blockFns.map((f) => f(INW));
  const inner = blocks.reduce((a, b) => a + b.h, 0);
  const H = PAD + headH + 0.35 + inner + PAD - 0.25;
  const x = colX[col], y = cursor[col];
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w: COLW, h: H, rectRadius: 0.3, fill: { color: C.background1 },
    line: { color: C.accent4, width: 2.5 }, objectName: `Section ${num} frame` });
  s.addShape(pres.shapes.OVAL, { x: x + PAD, y: y + PAD, w: headH, h: headH, fill: { color: C.accent1 }, line: { type: "none" }, objectName: `Section ${num} badge` });
  s.addText(String(num), { x: x + PAD, y: y + PAD, w: headH, h: headH, fontSize: 42, bold: true, color: C.background1, align: "center", valign: "middle", margin: 0, isTextBox: true });
  s.addText(title, { x: x + PAD + headH + 0.35, y: y + PAD, w: INW - headH - 0.35, h: headH, fontSize: HEAD, bold: true, color: C.text1, valign: "middle", margin: 0, isTextBox: true, objectName: `Section ${num} title` });
  let yy = y + PAD + headH + 0.35;
  blocks.forEach((b) => { b.draw(x + PAD, yy); yy += b.h; });
  cursor[col] = y + H + 0.5;
  console.log(`section ${num}: h=${H.toFixed(2)}  bottom=${(y + H).toFixed(2)}`);
}

// ---------------------------------------------------------------- header
s.addShape(pres.shapes.RECTANGLE, { x: 0, y: 0, w: PW, h: HDR_H, fill: { color: C.text1 }, line: { type: "none" }, objectName: "Header band" });
const LOGO = 4.4;
[["XENON logo", M], ["ISEE / Nagoya U. logo", PW - M - LOGO]].forEach(([lab, x]) => {
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: (HDR_H - LOGO) / 2, w: LOGO, h: LOGO, rectRadius: 0.2,
    fill: { color: C.background1, transparency: 88 }, line: { color: C.background1, width: 1.5, dashType: "dash" }, objectName: lab + " placeholder" });
  s.addText(lab + "\n(insert)", { x, y: (HDR_H - LOGO) / 2, w: LOGO, h: LOGO, fontSize: 24, color: C.background1, align: "center", valign: "middle", isTextBox: true });
});
const TX = M + LOGO + 0.6, TW = PW - 2 * TX;
s.addText("Measuring the Radon-222 Background in the XENONnT Dark Matter Detector", {
  x: TX, y: 0.45, w: TW, h: 3.1, fontSize: 86, bold: true, color: C.background1, align: "center", valign: "middle",
  margin: 0, isTextBox: true, objectName: "Title" });
s.addText([
  t("Jie Tang", { bold: true }), t(", on behalf of the XENON Collaboration", { breakLine: true }),
  t("Institute for Space-Earth Environmental Research (ISEE), Nagoya University"),
], { x: TX, y: 3.65, w: TW, h: 1.75, fontSize: 38, color: C.background1, align: "center", valign: "top", margin: 0,
  isTextBox: true, paraSpaceAfter: 6, objectName: "Authors" });

// ---------------------------------------------------------------- footer
s.addShape(pres.shapes.RECTANGLE, { x: 0, y: PH - FTR_H, w: PW, h: FTR_H, fill: { color: C.text1 }, line: { type: "none" }, objectName: "Footer band" });
s.addText([
  t("KMI2026", { bold: true }), t("  |  7th KMI International Symposium, Quest for the Origin of Particles and the Universe  |  Nagoya University, 9-11 November 2026", { breakLine: true }),
  t("[1] XENON Collaboration, Phys. Rev. X 15, 031079 (2025)    [2] XENON Collaboration, Phys. Rev. Lett. 135, 221003 (2025)    SR0/SR1 figures: Y.-T. Lin (XENON)", { fontSize: 20 }),
], { x: M, y: PH - FTR_H, w: PW * 0.72, h: FTR_H, fontSize: 25, color: C.background1, valign: "middle", margin: 0, isTextBox: true, paraSpaceAfter: 4, objectName: "Footer text" });
s.addText([t("Jie Tang", { bold: true, breakLine: true }), t("davidtang@nagoya-u.jp")], {
  x: PW - M - 7, y: PH - FTR_H, w: 7, h: FTR_H, fontSize: 25, color: C.background1, align: "right", valign: "middle", margin: 0, isTextBox: true, objectName: "Contact" });

const ORANGE = C.accent2, TEAL = C.accent1;
const credit = " Figure: Y.-T. Lin (XENON).";

// ================================================================ LEFT COLUMN
section(0, 1, "Radon in XENONnT", [
  para("XENONnT searches for dark matter with 5.9 t of liquid xenon in a dual-phase time projection chamber (TPC), deep underground at LNGS, Italy."),
  row([
    [0.42, [image("tpc_compact.png", { alt: "Dual-phase xenon TPC schematic" })]],
    [0.58, [bullets([
      [t("S1: ", { bold: true }), t("prompt scintillation light from the interaction")],
      [t("S2: ", { bold: true }), t("electrons drift up, are extracted into the gas and emit delayed light")],
      [t("Position: ", { bold: true }), t("x, y from the S2 light pattern; depth from the S1 to S2 delay")],
      [t("AFT: ", { bold: true }), t("fraction of S1 photoelectrons (PE) on the top PMTs; indicates depth from S1 alone")],
    ], { after: 0 })]],
  ], { vcenter: true, after: 0.35 }),
  sub("Why radon-222?"),
  bullets([
    "Continuously emanated from detector materials, it mixes into the xenon",
    [t("Its daughter "), sup("214"), t("Pb", { bold: true }), t(" beta decays inside the target: a "), t("leading low-energy background", { bold: true }), t(" for dark matter searches")],
    [t("The chain also emits "), t("mono-energetic alphas", { bold: true }), t(" at 5.6 to 7.8 MeV; counting them measures the radon level in situ")],
  ]),
  image("decay_chain_h.png", { alt: "Radon-222 and radon-220 decay chains" }),
  caption("Energies are decay Q-values. The thoron chain gives independent peaks as a stable reference."),
]);

section(0, 2, "Two cross-checked reconstruction modes", [
  row([
    [0.5, [sub("S1 + S2 (standard)", ORANGE), bullets([
      "Precise 3D position (S2 + drift time)",
      "S2 removes gammas, mis-paired and multiple-scatter events",
      [t("BiPo ("), sup("214"), t("Po) events not usable")],
    ], { after: 0 })]],
    [0.5, [sub("S1-only", TEAL), bullets([
      "Depth via AFT, radius via S1 light pattern (neural network)",
      [t("Keeps BiPo events: measures "), sup("214"), t("Po")],
      "Works without S2 (e.g. anode off)",
      "Volume calibrated with S1 + S2",
    ], { after: 0 })]],
  ]),
]);

section(0, 3, "Light-collection correction", [
  para([t("The "), t("light collection efficiency (LCE)", { bold: true }), t(" depends on where the alpha decays. Empirical corrections from "), sup("222"), t("Rn calibration data turn each curved band into a sharp line.")]),
  row([
    [0.5, [sub("S1 + S2: two steps", ORANGE), bullets([
      [t("Depth: ", { bold: true }), t("S1 vs AFT")],
      [t("Radius: ", { bold: true }), t("S1 vs r"), sup("2"), t(" from S2")],
    ], { after: 0 })]],
    [0.5, [sub("S1-only: four steps", TEAL), bullets([
      [t("Rough depth", { bold: true }), t(" to select "), sup("210"), t("Po")],
      [t("Radius: ", { bold: true }), sup("210"), t("Po on the walls fixes r"), sup("2")],
      [t("S1 vs corrected r"), sup("2")],
      [t("Final depth: ", { bold: true }), t("S1 vs AFT")],
    ], { after: 0 })]],
  ], { after: 0.1 }),
  row([
    [0.5, [image("B1_s1s2_corr.png", { alt: "S1+S2 correction before and after", name: "B1 S1+S2 correction" })]],
    [0.5, [image("B2_s1only_corr.png", { alt: "S1-only correction before and after", name: "B2 S1-only correction" })]],
  ], { after: 0.1 }),
  caption("AFT vs S1 for each mode, before (left) and after (right) correction, ²²²Rn calibration, SR1." + credit),
]);

// ================================================================ RIGHT COLUMN
section(1, 4, "Peak fit and rate", [
  row([
    [0.5, [image("B4_fits.png", { alt: "Multi-Gaussian fit of alpha peaks", name: "B4 fits" }),
      caption("One day of data. Top: S1 + S2; bottom: S1-only." + credit)]],
    [0.5, [bullets([
      "One Gaussian per alpha peak in the corrected S1 spectrum",
      [t("Peaks: "), sup("222"), t("Rn, "), sup("218"), t("Po, "), sup("220"), t("Rn, "), sup("216"), t("Po; S1-only adds "), sup("214"), t("Po")],
      [sup("212"), t("Bi alpha is negligible and not fitted")],
      [t("Rate", { bold: true }), t(" = counts within ±3σ / (acceptance × live time × fiducial mass), in μBq/kg")],
      "One rate per day of data",
    ], { after: 0 })]],
  ], { vcenter: true }),
]);

section(1, 5, "Results from SR0 and SR1", [
  row([
    [0.56, [image("PRX_fig6.png", { alt: "Radon-222 activity in XENONnT over time", name: "PRX Fig 6" }), caption("²²²Rn activity in XENONnT over time, from [1].")]],
    [0.44, [sub("Radon level in XENONnT"), bullets([
      [sup("222"), t("Rn: 3.6 → 1.9 μBq/kg (SR0) → "), t("0.90 μBq/kg", { bold: true }), t(" with cryogenic distillation [1]")],
      "About 5 times lower than other multi-tonne xenon detectors [1]",
      [sup("214"), t("Pb background reduced by 1.9 (1.7) in SR1a (SR1b) vs SR0 [2]")],
    ], { after: 0 })]],
  ], { vcenter: true, after: 0.25 }),
  row([
    [0.5, [image("B6_rateratio.png", { alt: "Radon-222 rate and low-energy rate", name: "B6 rate ratio" }),
      caption("²²²Rn calibration, SR1." + credit)]],
    [0.5, [sub("Alpha rate predicts the background"), bullets([
      [t("Fit of the rates: "), t("[10, 70] keV = a × ", { bold: true }), sup("222"), t("Rn + b", { bold: true })],
      [t("a", { bold: true }), t(": low-energy "), sup("214"), t("Pb events per radon decay; "), t("b", { bold: true }), t(": radon-independent background")],
      [t("The radon spike separates a from b: the ratio drops to a while radon dominates, then recovers as radon decays")],
    ], { after: 0 })]],
  ], { vcenter: true, after: 0.25 }),
  row([
    [0.5, [image("B5_ratios.png", { alt: "Isotope ratio comparison", name: "B5 ratios" }),
      caption("Change of each ratio relative to the ²²²Rn calibration." + credit)]],
    [0.5, [sub("The two modes agree"), bullets([
      [sup("218"), t("Po/"), sup("222"), t("Rn agrees within ~1% between "), t("S1 + S2 and S1-only", { bold: true }), t(" in every run")],
    ], { after: 0 })]],
  ], { vcenter: true }),
]);

section(1, 6, "Summary and outlook", [
  bullets([
    [t("Alpha spectroscopy with two cross-checked modes measures radon-222 "), t("in situ", { italic: true }), t(" and predicts the "), sup("214"), t("Pb background")],
  ], { after: 0.2 }),
  sub("Next: SR3"),
  bullets([
    "New detector configuration with an improved electric field",
    [t("Expected "), t("longer electron lifetime", { bold: true }), t(" and "), t("better energy resolution", { bold: true }), t(", leading to "), t("better sensitivity", { bold: true })],
    [t("Our radon framework is already applied to SR3 commissioning. "), t("Stay tuned!", { bold: true, color: C.accent1 })],
  ], { after: 0 }),
]);

console.log("bottoms:", cursor.map((c) => (c - 0.5).toFixed(2)), "limit:", (PH - FTR_H - 0.5).toFixed(2));

(async () => {
  const out = path.join(__dirname, "Rn222_XENONnT_KMI2026_poster.pptx");
  await pres.writeFile({ fileName: out });
  execSync(`python3 ${path.join(__dirname, "fix_ppr.py")} ${out}`);
  await applyTheme(out, THEME);
  console.log("wrote", out);
})();
