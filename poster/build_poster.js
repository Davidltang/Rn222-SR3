// A0 portrait poster: Measuring the Radon-222 Background in the XENONnT Dark Matter Detector
const pptxgen = require("pptxgenjs");
const path = require("path");
const { applyTheme } = require("/root/.claude/skills/synced/f6ad0e85-3217-498c-aa06-38e6c8806359_0ff8945c-f50f-42d2-854f-a9da0a03ddaa/pptx/scripts/apply_theme.js");

const IMG = (f) => path.join(__dirname, "img", f);
const imgSize = {
  "decay_chain.png": [3547, 1688], "tpc_schematic.png": [2851, 2023],
  "b1.png": [800, 578], "b2.png": [800, 289], "b3.png": [800, 307],
  "b4.png": [600, 596], "b5.png": [600, 269], "b6.png": [600, 537],
  "decay_chain_poster.png": [2643, 1902], "tpc_schematic_poster.png": [2679, 1931],
};

const THEME = {
  name: "Rn222 Poster",
  headFontFace: "Calibri",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "1B2A4A", lt1: "FFFFFF", dk2: "1B2A4A", lt2: "F2F5F8",
    accent1: "1F8A8A", accent2: "D9631E", accent3: "7A8291", accent4: "3D78B8",
    accent5: "E8B23A", accent6: "D5ECEC", hlink: "1F8A8A", folHlink: "7A8291",
  },
};
const HEX = { navy: "1B2A4A", teal: "1F8A8A", orange: "D9631E", grey: "5E6675", panel: "F2F5F8", tealL: "D5ECEC", white: "FFFFFF" };

const pres = new pptxgen();
pres.defineLayout({ name: "A0P", width: 33.11, height: 46.81 });
pres.layout = "A0P";
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "Measuring the Radon-222 Background in the XENONnT Dark Matter Detector";
pres.author = "Jie Tang";
const C = pres.SchemeColor;

pres.addSection({ title: "Poster" });
const s = pres.addSlide({ sectionTitle: "Poster" });
s.background = { color: HEX.white };

// ---------------------------------------------------------------- geometry
const PW = 33.11, PH = 46.81, M = 0.9, GAP = 0.6;
const COLW = (PW - 2 * M - 2 * GAP) / 3;
const colX = [0, 1, 2].map((i) => M + i * (COLW + GAP));
const PAD = 0.4;               // inside a section panel
const INW = COLW - 2 * PAD;    // inner width
const BODY = 25, CAP = 19, SUB = 28, HEAD = 44;
const PT = 1 / 72;

// rough Calibri line estimate: average glyph width ~0.48 em
function nLines(text, size, width) {
  const cpl = Math.max(10, Math.floor(width / (size * PT * 0.45)));
  return text.split("\n").reduce((n, para) => n + Math.max(1, Math.ceil(para.length / cpl)), 0);
}
const lineH = (size, sp = 1.18) => size * PT * sp;

// ---------------------------------------------------------------- header
const HDR_H = 6.2;
s.addShape(pres.shapes.RECTANGLE, { x: 0, y: 0, w: PW, h: HDR_H, fill: { color: C.text1 }, line: { type: "none" }, objectName: "Header background" });
s.addText("Measuring the Radon-222 Background in the XENONnT Dark Matter Detector", {
  x: M, y: 0.75, w: PW - 2 * M - 6.6, h: 3.1, fontSize: 92, bold: true, color: C.background1,
  valign: "top", margin: 0, isTextBox: true, fontFace: "Calibri", objectName: "Title",
});
s.addText([
  { text: "Jie Tang", options: { bold: true, breakLine: true } },
  { text: "Institute for Space-Earth Environmental Research (ISEE), Nagoya University", options: { breakLine: true } },
  { text: "on behalf of the XENON Collaboration", options: { italic: true } },
], { x: M, y: 4.0, w: PW - 2 * M - 6.6, h: 1.8, fontSize: 36, color: C.background1, valign: "top", margin: 0, isTextBox: true, paraSpaceAfter: 4, objectName: "Authors" });
// logo placeholders
[["XENON logo", 0.75], ["ISEE / Nagoya U. logo", 3.35]].forEach(([lab, y]) => {
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: PW - M - 5.8, y, w: 5.8, h: 2.2, rectRadius: 0.15,
    fill: { color: C.background1, transparency: 85 }, line: { color: C.background1, width: 1.5, dashType: "dash" }, objectName: lab + " placeholder" });
  s.addText(lab + " (insert)", { x: PW - M - 5.8, y, w: 5.8, h: 2.2, fontSize: 24, color: C.background1, align: "center", valign: "middle", isTextBox: true });
});

// ---------------------------------------------------------------- flow layout
const cursor = [HDR_H + 0.8, HDR_H + 0.8, HDR_H + 0.8];

// a section = list of blocks; each block returns {h, draw(x,y)}
function para(runs, opts = {}) {
  const size = opts.size || BODY;
  const width = opts.w || INW;
  const plain = (Array.isArray(runs) ? runs.map((r) => r.text + (r.options && r.options.breakLine ? "\n" : "")).join("") : runs);
  const nl = nLines(plain.replace(/\n$/, ""), size, width - 0.1);
  const h = nl * lineH(size) + 0.12 + (opts.extra || 0);
  return { h: h + (opts.after ?? 0.25), draw: (x, y) => s.addText(runs, {
    x, y, w: width, h, fontSize: size, color: opts.color || C.text1, margin: 0, valign: "top",
    isTextBox: true, align: "left", lineSpacingMultiple: 1.08, italic: opts.italic, bold: opts.bold, paraSpaceAfter: opts.psa || 0,
  }) };
}
function bullets(items, opts = {}) {
  const size = opts.size || BODY;
  const runs = [];
  let nl = 0;
  items.forEach((it, i) => {
    const parts = Array.isArray(it) ? it : [{ text: it }];
    parts.forEach((p, j) => {
      const o = Object.assign({}, p.options || {});
      if (j === 0) o.bullet = { indent: 26 };
      if (j === parts.length - 1 && i < items.length - 1) o.breakLine = true;
      o.paraSpaceAfter = 8;
      runs.push({ text: p.text, options: o });
    });
    nl += nLines(parts.map((p) => p.text).join(""), size, INW - 0.55);
  });
  const h = nl * lineH(size) + items.length * 8 * PT + 0.15;
  return { h: h + 0.25, draw: (x, y) => s.addText(runs, { x, y, w: INW, h, fontSize: size, color: C.text1, margin: 0, valign: "top", isTextBox: true, lineSpacingMultiple: 1.08 }) };
}
function sub(text, color) {
  const h = lineH(SUB) + 0.08;
  return { h: h + 0.12, draw: (x, y) => s.addText(text, { x, y, w: INW, h, fontSize: SUB, bold: true, color: color || C.accent1, margin: 0, valign: "top", isTextBox: true }) };
}
function image(file, opts = {}) {
  const [pw, ph] = imgSize[file];
  const w = opts.w || INW;
  const h = w * ph / pw;
  return { h: h + (opts.after ?? 0.12), draw: (x, y) => s.addImage({ path: IMG(file), x: x + (INW - w) / 2, y, w, h, altText: opts.alt || file, objectName: opts.name || file }) };
}
function caption(text) {
  const h = nLines(text, CAP, INW) * lineH(CAP, 1.15) + 0.08;
  return { h: h + 0.3, draw: (x, y) => s.addText(text, { x, y, w: INW, h, fontSize: CAP, italic: true, color: C.accent3, margin: 0, valign: "top", isTextBox: true }) };
}
function placeholder(label, hgt) {
  return { h: hgt + 0.12, draw: (x, y) => {
    s.addShape(pres.shapes.RECTANGLE, { x, y, w: INW, h: hgt, fill: { color: C.background1 }, line: { color: C.accent3, width: 2, dashType: "dash" }, objectName: label });
    s.addText(label, { x, y, w: INW, h: hgt, fontSize: 26, color: C.accent3, align: "center", valign: "middle", isTextBox: true });
  } };
}
function callout(text) {
  const size = 30;
  const h = nLines(text, size, INW - 0.8) * lineH(size) + 0.55;
  return { h: h + 0.3, draw: (x, y) => s.addText(text, { x, y, w: INW, h, fontSize: size, bold: true, color: C.background1,
    fill: { color: C.accent1 }, align: "center", valign: "middle", margin: [10, 18, 10, 18], isTextBox: true, shape: pres.shapes.ROUNDED_RECTANGLE, rectRadius: 0.15 }) };
}

function section(col, num, title, blocks) {
  const headH = 0.9;
  const inner = blocks.reduce((a, b) => a + b.h, 0);
  const H = PAD + headH + 0.3 + inner + PAD - 0.2;
  const x = colX[col], y = cursor[col];
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w: COLW, h: H, rectRadius: 0.25, fill: { color: C.background2 }, line: { type: "none" }, objectName: `Section ${num} panel` });
  s.addShape(pres.shapes.OVAL, { x: x + PAD, y: y + PAD, w: headH, h: headH, fill: { color: C.accent1 }, line: { type: "none" }, objectName: `Section ${num} badge` });
  s.addText(String(num), { x: x + PAD, y: y + PAD, w: headH, h: headH, fontSize: 38, bold: true, color: C.background1, align: "center", valign: "middle", margin: 0, isTextBox: true });
  s.addText(title, { x: x + PAD + headH + 0.3, y: y + PAD, w: INW - headH - 0.3, h: headH, fontSize: HEAD, bold: true, color: C.text1, valign: "middle", margin: 0, isTextBox: true, objectName: `Section ${num} title` });
  let yy = y + PAD + headH + 0.3;
  blocks.forEach((b) => { b.draw(x + PAD, yy); yy += b.h; });
  cursor[col] = y + H + 0.55;
  console.log('section', num, 'h', H.toFixed(2));
  return H;
}

const sup = (t) => ({ text: t, options: { superscript: true } });
const t = (txt, o = {}) => ({ text: txt, options: o });

// ================================================================ COLUMN 1
section(0, 1, "Why radon?", [
  para("XENONnT searches for dark matter with 5.9 t of liquid xenon (LXe) deep underground at LNGS. At this sensitivity, the leading backgrounds come from inside the detector itself."),
  bullets([
    [t("Radon-222 ", { bold: true }), t("is continuously emanated from detector materials and mixes into the xenon.")],
    [t("Its daughter "), sup("214"), t("Pb ", { bold: true }), t("beta decays with no accompanying signal, producing low-energy electron recoils in the dark matter search region.")],
    [t("The same chain emits "), t("mono-energetic alphas", { bold: true }), t(" at 5.6 to 7.8 MeV. Counting them measures the radon activity in situ.")],
  ]),
  image("decay_chain_poster.png", { w: INW * 0.86, alt: "Radon-222 and radon-220 decay chains" }),
  caption("Teal: alpha peaks used in the rate fit; orange: the ²¹⁴Pb background; grey: not fitted. Alpha energies are decay Q-values (alpha + nuclear recoil, both deposited in the xenon). The radon-220 (thoron) chain provides independent peaks that serve as a stable reference."),
]);

section(0, 2, "Detector and signals", [
  para([t("XENONnT is a dual-phase time projection chamber (TPC) viewed by photomultiplier (PMT) arrays at top and bottom. An interaction gives prompt scintillation light ("), t("S1", { bold: true }), t("); its ionization electrons drift up, are extracted into the gas and produce a delayed, amplified light signal ("), t("S2", { bold: true }), t(").")]),
  image("tpc_schematic_poster.png", { w: INW * 0.86, alt: "Schematic of a dual-phase xenon TPC" }),
  caption("Schematic, not to scale. The S2 pattern on the top array gives x and y, the S1 to S2 delay gives the depth, and the S1 area fraction top (AFT) gives the depth from S1 alone."),
]);

// ================================================================ COLUMN 2
section(0, 3, "Two reconstruction modes", [
  sub("S1 + S2 (standard)", C.accent2),
  para("Each alpha S1 is paired with its S2. Position comes from the S2 (x, y) and the drift time (z), independently of the S1 light being measured. The S2 also rejects gamma events, which have a much larger S2/S1 than alphas, unphysical S1-S2 pairings and multiple scatters, and defines a precise fiducial volume."),
  sub("S1-only", C.accent1),
  para([t("Only the S1 is used: depth from AFT, radius from a neural network trained on the S1 light pattern (S1 CNN). This keeps the BiPo decays ("), sup("214"), t("Bi beta followed by "), sup("214"), t("Po alpha within 164 μs), whose S2 is confused by the preceding beta, and works when no S2 is available, for example with the anode off. Its fiducial mass is calibrated by matching events to their S1 + S2 positions.")]),
]);

section(1, 4, "Light-collection correction", [
  para("The S1 light collected from an alpha depends on where it decays. Empirical corrections, derived from radon-222 calibration data, turn each curved alpha band into a narrow vertical line."),
  sub("S1 + S2: two steps", C.accent2),
  bullets([
    [t("Depth: ", { bold: true }), t("S1 area vs AFT (drift time gives the same result).")],
    [t("Radius: ", { bold: true }), t("residual dependence on r"), sup("2"), t(" from the S2 position.")],
  ]),
  image("b1.png", { w: INW * 0.9, alt: "S1+S2 alpha selection and correction", name: "B1 preview (replace with original)" }),
  caption("S1 + S2, radon-222 calibration, SR1. Top: S2 vs S1 and AFT vs S1 before correction. Bottom: after correction, each isotope is a separate line. Figure: Y.-T. Lin (XENON)."),
  sub("S1-only: four steps", C.accent1),
  bullets([
    [t("Preliminary depth correction ", { bold: true }), t("(area vs AFT) to select "), sup("210"), t("Po.")],
    [t("Radius: ", { bold: true }), sup("210"), t("Po sits on the walls (plate-out), giving a fixed reference to correct the S1 CNN r"), sup("2"), t(" vs AFT.")],
    [t("Area vs corrected r"), sup("2"), t(".", {})],
    [t("Final depth correction ", { bold: true }), t("(area vs AFT) in a tighter volume.")],
  ]),
  image("b2.png", { alt: "S1-only alpha correction", name: "B2 preview (replace with original)" }),
  caption("S1-only, radon-222 calibration, SR1: AFT vs S1 before (left) and after (right) correction. Figure: Y.-T. Lin (XENON)."),
]);

// ================================================================ COLUMN 3
section(1, 5, "Multi-Gaussian fit", [
  image("b4.png", { w: INW * 0.78, alt: "Multi-Gaussian fit of alpha peaks", name: "B4 preview (replace with original)" }),
  caption("Fit to one day of data. Top: S1 + S2; bottom: S1-only, which also recovers the ²¹⁴Po peak. Figure: Y.-T. Lin (XENON)."),
  bullets([
    "Data are grouped into 1-day chunks; peak widths are linked through the square root of the peak position.",
    "Each rate is the count within ±3σ of its peak, corrected for the Gaussian acceptance.",
    [t("The "), sup("212"), t("Bi alpha branch lies between the "), sup("218"), t("Po and "), sup("220"), t("Rn peaks; its contribution is negligible and is not fitted.")],
  ]),
]);

section(2, 6, "Results from SR0 and SR1", [
  sub("Alpha rate tracks the low-energy background"),
  image("b6.png", { w: INW * 0.9, alt: "Radon-222 rate and low-energy rate during calibration", name: "B6 preview (replace with original)" }),
  caption("Radon-222 calibration, SR1. The rate of [10, 70] keV events follows the radon-222 alpha rate; their ratio links the alpha measurement to the ²¹⁴Pb background. Figure: Y.-T. Lin (XENON)."),
  sub("Consistent isotope ratios"),
  image("b5.png", { alt: "Isotope ratio comparison", name: "B5 preview (replace with original)" }),
  caption("²¹⁸Po/²²²Rn and ²¹⁴Po/²²²Rn across SR0, SR1 and calibration, for both modes. The spread in ²¹⁴Po/²²²Rn reflects plate-out and is taken as a systematic. Figure: Y.-T. Lin (XENON)."),
  sub("Radon reduction in XENONnT"),
  placeholder("PRX 15, 031079 (2025), Fig. 6 (insert)", 5.2),
  caption("Published radon-222 activity: 3.6 μBq/kg without radon removal, 1.87 μBq/kg in SR0 and down to 0.90 μBq/kg with removal in both gas and liquid [1, 2]."),
]);

section(2, 7, "Outlook: SR3", [
  para("The same framework is now applied to SR3 commissioning data, taken in a new detector configuration including a revised electric field. The S1-only mode follows the radon level through periods without S2, with S1 + S2 cross-checks where both are available."),
]);

section(2, 8, "Summary", [
  callout("Count alphas at MeV energies to predict the background at keV energies"),
  bullets([
    [t("Radon-222 is measured "), t("in situ", { italic: true }), t(" through its mono-energetic alpha peaks, giving direct access to the "), sup("214"), t("Pb background.")],
    [t("Two complementary modes, S1 + S2 and S1-only, are cross-checked; S1-only adds "), sup("214"), t("Po and works without S2.")],
    "Demonstrated on SR0 and SR1; now extended to SR3.",
  ]),
]);

// ---------------------------------------------------------------- footer
const footY = PH - M - 1.25;
s.addText([
  t("References  ", { bold: true }),
  t("[1] XENON Collaboration, Phys. Rev. X 15, 031079 (2025).   [2] XENON Collaboration, Phys. Rev. Lett. 135, 221003 (2025).", { breakLine: true }),
  t("SR0/SR1 analysis and figures: Y.-T. Lin (MPIK) on behalf of the XENON Collaboration.   Contact: davidtang@nagoya-u.jp"),
], { x: M, y: footY, w: PW - 2 * M, h: 1.25, fontSize: 22, color: C.accent3, valign: "bottom", margin: 0, isTextBox: true, objectName: "Footer" });

console.log("column bottoms:", cursor.map((c) => c.toFixed(2)), "footer at", footY.toFixed(2));

(async () => {
  const out = path.join(__dirname, "Rn222_XENONnT_KMI2026_poster.pptx");
  await pres.writeFile({ fileName: out });
  require("child_process").execSync(`python3 ${path.join(__dirname, "fix_ppr.py")} ${out}`);
  await applyTheme(out, THEME);
  console.log("wrote", out);
})();
