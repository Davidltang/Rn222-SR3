import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle, Polygon
import numpy as np

plt.rcParams.update({
    "font.family": "Carlito",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Carlito",
    "mathtext.it": "Carlito:italic",
    "mathtext.bf": "Carlito:bold",
})

NAVY = "#1B2A4A"
TEAL = "#1F8A8A"
TEAL_L = "#D5ECEC"
ORANGE = "#D9631E"
ORANGE_L = "#FBE3D3"
GREY = "#7A8291"
GREY_L = "#ECEEF1"
GOLD = "#E8B23A"
OUT = "/home/claude/poster/fig/"


def iso(a, el):
    return r"$^{%s}$%s" % (a, el)


def box(ax, x, y, w, h, fc, ec, lw=2.0, r=0.08):
    p = FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                       boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw)
    ax.add_patch(p)


def arrow(ax, x0, y0, x1, y1, color=NAVY, lw=2.2, style="-|>", ms=18, ls="-"):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style, mutation_scale=ms,
                        color=color, lw=lw, linestyle=ls, shrinkA=0, shrinkB=0)
    ax.add_patch(a)


# ---------------------------------------------------------------- 1. decay chain
def decay_chain():
    fig, ax = plt.subplots(figsize=(15, 7.05))
    ax.set_xlim(0, 15)
    ax.set_ylim(-0.45, 6.6)
    ax.axis("off")
    W, H = 1.35, 0.95
    xs = [0.85, 3.4, 5.95, 8.5, 11.05, 13.6]

    def node(x, y, a, el, t, kind):
        fc, ec = {"alpha": (TEAL_L, TEAL), "bkg": (ORANGE_L, ORANGE),
                  "beta": (GREY_L, GREY), "end": ("white", GREY)}[kind]
        box(ax, x, y, W, H, fc, ec, lw=2.4 if kind in ("alpha", "bkg") else 1.6)
        ax.text(x, y + 0.13, iso(a, el), ha="center", va="center", fontsize=21,
                color=NAVY, fontweight="bold")
        ax.text(x, y - 0.25, t, ha="center", va="center", fontsize=13, color=GREY)

    def step(x0, x1, y, lab, kind):
        c = TEAL if kind == "a" else GREY
        arrow(ax, x0 + W / 2 + 0.04, y, x1 - W / 2 - 0.06, y, color=c)
        ax.text((x0 + x1) / 2, y + 0.24, lab, ha="center", va="bottom",
                fontsize=13, color=c, fontweight="bold" if kind == "a" else "normal")

    # radon-222 chain
    y1 = 5.0
    ax.text(0.2, 6.25, r"Radon-222 chain (uranium series)", fontsize=19, color=NAVY,
            fontweight="bold", va="center")
    node(xs[0], y1, 222, "Rn", "3.82 d", "alpha")
    node(xs[1], y1, 218, "Po", "3.10 min", "alpha")
    node(xs[2], y1, 214, "Pb", "26.8 min", "bkg")
    node(xs[3], y1, 214, "Bi", "19.9 min", "beta")
    node(xs[4], y1, 214, "Po", r"164 $\mu$s", "alpha")
    node(xs[5], y1, 210, "Pb", "22.2 yr", "end")
    step(xs[0], xs[1], y1, r"$\alpha$ 5.59 MeV", "a")
    step(xs[1], xs[2], y1, r"$\alpha$ 6.11 MeV", "a")
    step(xs[2], xs[3], y1, r"$\beta$", "b")
    step(xs[3], xs[4], y1, r"$\beta$", "b")
    step(xs[4], xs[5], y1, r"$\alpha$ 7.83 MeV", "a")
    ax.text(xs[2], y1 - 0.72, "low-energy background\nin the dark matter search",
            ha="center", va="top", fontsize=13, color=ORANGE, fontweight="bold")
    # BiPo bracket
    bx0, bx1, by = xs[3] - W / 2, xs[4] + W / 2, y1 - 0.62
    ax.plot([bx0, bx0, bx1, bx1], [by + 0.08, by, by, by + 0.08], color=NAVY, lw=1.4)
    ax.text((bx0 + bx1) / 2, by - 0.1, "\"BiPo\": two decays within\n" + r"~0.2 ms of each other",
            ha="center", va="top", fontsize=13, color=NAVY)

    # thoron chain
    y2 = 2.05
    ax.text(0.2, 3.3, r"Radon-220 chain (thorium series, stable reference)", fontsize=19,
            color=NAVY, fontweight="bold", va="center")
    node(xs[0], y2, 220, "Rn", "55.6 s", "alpha")
    node(xs[1], y2, 216, "Po", "0.145 s", "alpha")
    node(xs[2], y2, 212, "Pb", "10.6 h", "beta")
    node(xs[3], y2, 212, "Bi", "60.6 min", "beta")
    step(xs[0], xs[1], y2, r"$\alpha$ 6.40 MeV", "a")
    step(xs[1], xs[2], y2, r"$\alpha$ 6.91 MeV", "a")
    step(xs[2], xs[3], y2, r"$\beta$", "b")
    # Bi212 branches
    xb = xs[4] - 0.2
    arrow(ax, xs[3] + W / 2 + 0.04, y2 + 0.15, xb - 0.95, y2 + 0.55, color=GREY)
    arrow(ax, xs[3] + W / 2 + 0.04, y2 - 0.15, xb - 0.95, y2 - 0.55, color=GREY)
    ax.text(xb - 0.85, y2 + 0.6, r"$\alpha$ (36%) $\rightarrow$ $^{208}$Tl",
            ha="left", va="center", fontsize=13, color=GREY)
    ax.text(xb - 0.85, y2 - 0.6, r"$\beta$ (64%) $\rightarrow$ $^{212}$Po ($\alpha$ 8.95 MeV, 0.3 $\mu$s)",
            ha="left", va="center", fontsize=13, color=GREY)

    # legend
    ly = 0.45
    for i, (fc, ec, lab) in enumerate([(TEAL_L, TEAL, "alpha emitter (measured peak)"),
                                       (ORANGE_L, ORANGE, "beta background of interest"),
                                       (GREY_L, GREY, "not used in the rate fit")]):
        x = [0.35, 4.5, 8.55][i]
        box(ax, x + 0.25, ly, 0.5, 0.36, fc, ec, lw=1.8, r=0.05)
        ax.text(x + 0.65, ly, lab, va="center", fontsize=14, color=NAVY)
    ax.text(0.35, -0.15, "Alpha energies are decay Q-values (alpha + nuclear recoil, both deposited in the liquid xenon).",
            ha="left", va="center", fontsize=13, color=GREY)
    fig.savefig(OUT + "decay_chain.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------- 2. TPC schematic
def tpc():
    fig, ax = plt.subplots(figsize=(11.6, 8.5))
    ax.set_xlim(0, 11.6)
    ax.set_ylim(0, 8.5)
    ax.set_aspect("equal")
    ax.axis("off")
    x0, x1 = 0.8, 5.4
    yb, yt = 0.7, 7.9
    y_cath, y_gate, y_liq, y_anode = 1.25, 6.35, 6.6, 6.85
    # vessel
    ax.add_patch(Rectangle((x0, yb), x1 - x0, yt - yb, fc="white", ec=NAVY, lw=2.5))
    ax.add_patch(Rectangle((x0, y_cath - 0.05), x1 - x0, y_liq - y_cath + 0.05,
                           fc="#DCE9F7", ec="none"))
    ax.add_patch(Rectangle((x0, y_liq), x1 - x0, yt - 0.55 - y_liq, fc="#F4F8FC", ec="none"))
    ax.plot([x0, x1], [y_liq, y_liq], color="#3D78B8", lw=2)
    # PMT arrays
    for yy in (yt - 0.3, yb + 0.3):
        for xc in np.linspace(x0 + 0.35, x1 - 0.35, 8):
            ax.add_patch(Circle((xc, yy), 0.22, fc=GOLD, ec=NAVY, lw=1.2))
    # electrodes
    for yy, lab in ((y_anode, "anode"), (y_gate, "gate"), (y_cath, "cathode")):
        ax.plot([x0, x1], [yy, yy], color=GREY, lw=1.6, ls=(0, (2, 2)))
        ax.text(x0 - 0.12, yy, lab, ha="right", va="center", fontsize=13, color=GREY)
    ax.text(x1 + 0.12, y_anode + 0.25, "gas xenon", fontsize=13, color=NAVY, va="center")
    ax.text(x1 + 0.12, 3.6, "liquid xenon", fontsize=13, color=NAVY, va="center")
    ax.text(x1 + 0.12, yt - 0.3, "top PMT array", fontsize=13, color=NAVY, va="center")
    ax.text(x1 + 0.12, yb + 0.3, "bottom PMT array", fontsize=13, color=NAVY, va="center")
    # extraction field bracket
    ex = x1 + 0.12
    ax.plot([ex, ex + 0.1, ex + 0.1, ex], [y_gate, y_gate, y_anode, y_anode], color=ORANGE, lw=1.5)
    ax.text(ex + 0.22, (y_gate + y_anode) / 2 - 0.02, "extraction field", fontsize=12.5, color=ORANGE, va="center")
    # drift field arrow
    arrow(ax, x0 + 0.3, y_gate - 0.3, x0 + 0.3, y_cath + 0.3, color=GREY, lw=1.5, ms=13)
    ax.text(x0 + 0.42, 3.2, "drift field", rotation=90, fontsize=12, color=GREY, va="center")
    # interaction
    ix, iy = 3.4, 2.6
    star = Polygon([(ix + 0.32 * np.cos(t) * (1 if k % 2 == 0 else 0.45),
                     iy + 0.32 * np.sin(t) * (1 if k % 2 == 0 else 0.45))
                    for k, t in enumerate(np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, 11)[:-1])],
                   fc=ORANGE, ec=NAVY, lw=1)
    ax.add_patch(star)
    # S1 light rays
    for tx in (1.4, 2.4, 4.9):
        arrow(ax, ix, iy + 0.3, tx, yt - 0.55, color=TEAL, lw=1.4, ms=12)
    for bx in (2.3, 3.9, 4.9):
        arrow(ax, ix, iy - 0.3, bx, yb + 0.55, color=TEAL, lw=1.4, ms=12)
    ax.text(ix + 0.38, iy - 0.05, "S1", fontsize=20, color=TEAL, fontweight="bold", va="center")
    # electrons drifting
    arrow(ax, ix, iy + 0.35, ix, y_gate - 0.05, color="#3D78B8", lw=2, ms=15, ls=(0, (4, 3)))
    ax.text(ix + 0.12, 4.4, r"electrons drift up", rotation=90, fontsize=12, color="#3D78B8",
            ha="left", va="center")
    # S2 flash
    ax.add_patch(Circle((ix, (y_liq + y_anode) / 2 + 0.05), 0.2, fc="#FFE9A8", ec=ORANGE, lw=1.5))
    ax.text(ix + 0.3, (y_liq + y_anode) / 2 + 0.38, "S2", fontsize=20, color=ORANGE,
            fontweight="bold", va="center")

    # right panel: waveform
    wx0, wx1, wy0 = 7.6, 11.3, 5.2
    ax.text(wx0, 7.75, "Recorded signals", fontsize=17, color=NAVY, fontweight="bold")
    t = np.linspace(0, 1, 600)
    s1 = 0.55 * np.exp(-0.5 * ((t - 0.12) / 0.008) ** 2)
    s2 = 1.6 * np.exp(-0.5 * ((t - 0.78) / 0.04) ** 2)
    ax.plot(wx0 + t * (wx1 - wx0), wy0 + s1 + s2, color=NAVY, lw=1.8)
    ax.plot([wx0, wx1], [wy0, wy0], color=GREY, lw=1)
    ax.text(wx0 + 0.12 * (wx1 - wx0), wy0 + 0.72, "S1", ha="center", fontsize=15, color=TEAL, fontweight="bold")
    ax.text(wx0 + 0.78 * (wx1 - wx0) + 0.25, wy0 + 1.55, "S2", ha="left", fontsize=15, color=ORANGE, fontweight="bold")
    arrow(ax, wx0 + 0.12 * (wx1 - wx0), wy0 - 0.25, wx0 + 0.78 * (wx1 - wx0), wy0 - 0.25,
          color=GREY, lw=1.3, style="<|-|>", ms=10)
    ax.text(wx0 + 0.45 * (wx1 - wx0), wy0 - 0.55, "drift time = depth", ha="center", fontsize=13, color=GREY)
    ax.text(wx1, wy0 - 0.25, "time", ha="right", va="top", fontsize=11, color=GREY)

    # definitions
    ty = 3.95
    lines = [
        (TEAL, "S1", "prompt scintillation light"),
        (ORANGE, "S2", "delayed light from extracted electrons"),
    ]
    for i, (c, k, v) in enumerate(lines):
        ax.text(wx0, ty - i * 0.5, k, fontsize=15, color=c, fontweight="bold", va="center")
        ax.text(wx0 + 0.5, ty - i * 0.5, v, fontsize=14, color=NAVY, va="center")
    box(ax, (wx0 + wx1) / 2 - 0.05, 1.65, wx1 - wx0 + 0.3, 2.1, GREY_L, NAVY, lw=1.4, r=0.1)
    ax.text(wx0 - 0.05, 2.35, "AFT (area fraction top)", fontsize=15, color=NAVY, fontweight="bold", va="center")
    ax.text(wx0 - 0.05, 1.8, r"= S1 light on top array / total S1", fontsize=14, color=NAVY, va="center")
    ax.text(wx0 - 0.05, 1.3, "events near the top give a higher AFT,\nso AFT tracks depth using S1 alone",
            fontsize=13, color=GREY, va="center", linespacing=1.3)
    fig.savefig(OUT + "tpc_schematic.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------- 3. mode comparison
def modes():
    fig, ax = plt.subplots(figsize=(13, 6.2))
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 6.2)
    ax.axis("off")
    cols = [
        (3.05, ORANGE, ORANGE_L, "S1 + S2", "standard reconstruction", [
            ("Signals", "light (S1) and charge (S2)"),
            ("Position", "full 3D: x, y from S2 pattern,\nz from drift time"),
            ("Strength", "precise fiducial volume,\nstrong background rejection"),
            ("Limitation", r"BiPo ($^{214}$Po) events overlap,"+"\nnot measured"),
        ]),
        (9.95, TEAL, TEAL_L, "S1-only", "scintillation light only", [
            ("Signals", "light (S1) only"),
            ("Position", "depth from AFT,\nradius from S1 light pattern"),
            ("Strength", r"measures $^{214}$Po (BiPo);"+"\nworks when S2 is unavailable"),
            ("Limitation", "coarser position; volume\ncalibrated against S1 + S2"),
        ]),
    ]
    for cx, c, cl, title, sub, rows in cols:
        box(ax, cx, 3.05, 5.6, 5.6, "white", c, lw=2.6, r=0.18)
        ax.add_patch(FancyBboxPatch((cx - 2.8, 4.95), 5.6, 0.9,
                                    boxstyle="round,pad=0,rounding_size=0.18", fc=c, ec=c))
        ax.add_patch(Rectangle((cx - 2.8, 4.95), 5.6, 0.3, fc=c, ec=c))
        ax.text(cx, 5.5, title, ha="center", va="center", fontsize=24, color="white", fontweight="bold")
        ax.text(cx, 5.12, sub, ha="center", va="center", fontsize=14, color="white")
        for i, (k, v) in enumerate(rows):
            y = 4.4 - i * 1.05
            ax.text(cx - 2.55, y, k, fontsize=15, color=c, fontweight="bold", va="top")
            ax.text(cx - 1.0, y, v, fontsize=15, color=NAVY, va="top", linespacing=1.25)
    # cross-check link
    arrow(ax, 5.95, 3.05, 7.05, 3.05, color=NAVY, style="<|-|>", ms=16, lw=2.2)
    ax.text(6.5, 3.45, "cross-\ncheck", ha="center", va="bottom", fontsize=14, color=NAVY,
            fontweight="bold", linespacing=1.1)
    fig.savefig(OUT + "modes_comparison.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------- 4. correction flow
def correction():
    fig, ax = plt.subplots(figsize=(14, 5.2))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 5.2)
    ax.axis("off")
    steps = [
        ("1", "Depth", "S1 area vs AFT",
         "light collection changes\nwith depth: flatten the band", "curve"),
        ("2", "Decouple", r"r$^2$ vs AFT",
         r"remove the depth bias of the"+"\nreconstructed radius", "tilt"),
        ("3", "Radial", r"S1 area vs corrected r$^2$",
         "remove the remaining\nradial dependence", "radial"),
    ]
    cx0, dx = 2.15, 4.6
    for i, (n, name, axes_lab, desc, kind) in enumerate(steps):
        cx = cx0 + i * dx
        box(ax, cx, 2.6, 3.9, 4.9, "white", NAVY, lw=1.8, r=0.15)
        ax.add_patch(Circle((cx - 1.55, 4.55), 0.3, fc=TEAL, ec="none"))
        ax.text(cx - 1.55, 4.55, n, ha="center", va="center", fontsize=19, color="white", fontweight="bold")
        ax.text(cx - 1.1, 4.55, name, va="center", fontsize=21, color=NAVY, fontweight="bold")
        ax.text(cx, 3.95, axes_lab, ha="center", va="center", fontsize=14, color=GREY)
        # mini before/after sketch
        for j, (ox, lab) in enumerate(((cx - 0.95, "before"), (cx + 0.95, "after"))):
            ax.add_patch(Rectangle((ox - 0.7, 1.95), 1.4, 1.55, fc=GREY_L, ec=GREY, lw=1))
            yy = np.linspace(2.05, 3.4, 50)
            if j == 0:
                if kind == "curve":
                    xx = ox - 0.25 + 0.55 * ((yy - 2.05) / 1.35) ** 1.7
                elif kind == "tilt":
                    xx = ox - 0.4 + 0.75 * (yy - 2.05) / 1.35
                else:
                    xx = ox + 0.3 - 0.5 * ((yy - 2.05) / 1.35) ** 2.2
            else:
                xx = np.full_like(yy, ox)
            ax.fill_betweenx(yy, xx - 0.1, xx + 0.1, color=TEAL, alpha=0.75, lw=0)
            ax.text(ox, 1.78, lab, ha="center", va="top", fontsize=12, color=GREY)
        arrow(ax, cx - 0.2, 2.72, cx + 0.2, 2.72, color=NAVY, lw=1.6, ms=12)
        ax.text(cx, 1.05, desc, ha="center", va="center", fontsize=14, color=NAVY, linespacing=1.25)
        if i < 2:
            arrow(ax, cx + 2.0, 2.6, cx + dx - 2.0, 2.6, color=TEAL, lw=3, ms=22)
    fig.savefig(OUT + "correction_flow.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    import os
    os.makedirs(OUT, exist_ok=True)
    decay_chain(); tpc(); modes(); correction()


# ======================================================================
# Poster-scale versions: figure size = printed size on A0, fonts in true pt
# ======================================================================
def decay_chain_poster(width=9.2):
    Hh = 6.9
    fig = plt.figure(figsize=(width, Hh))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(0, Hh)
    ax.set_aspect("equal")
    ax.axis("off")
    BW, BH = 1.1, 0.72
    FS_ISO, FS_T, FS_L = 21, 15, 16

    def node(x, y, a, el, tt, kind):
        fc, ec = {"alpha": (TEAL_L, TEAL), "bkg": (ORANGE_L, ORANGE),
                  "beta": (GREY_L, GREY), "end": ("white", GREY)}[kind]
        box(ax, x, y, BW, BH, fc, ec, lw=2.2 if kind in ("alpha", "bkg") else 1.4, r=0.06)
        ax.text(x, y + 0.11, iso(a, el), ha="center", va="center", fontsize=FS_ISO, color=NAVY, fontweight="bold")
        ax.text(x, y - 0.2, tt, ha="center", va="center", fontsize=FS_T - 2, color=GREY)

    def adown(x, y0, y1, lab, measured=True):
        c = TEAL if measured else GREY
        arrow(ax, x, y0 - BH / 2 - 0.04, x, y1 + BH / 2 + 0.05, color=c, lw=2.2, ms=16)
        ax.text(x + 0.1, (y0 + y1) / 2, lab, ha="left", va="center", rotation=90, fontsize=FS_L, color=c,
                fontweight="bold" if measured else "normal")

    def bup(x0, y0, x1, y1, lab="β", dy=0.1):
        arrow(ax, x0 + BW / 2 - 0.05, y0 + dy, x1 - BW / 2 + 0.05, y1 - dy, color=GREY, lw=1.8, ms=14)
        ax.text((x0 + x1) / 2 - 0.1, (y0 + y1) / 2 + 0.12, lab, ha="right", va="bottom", fontsize=FS_L, color=GREY)

    r0, r1, r2 = 5.55, 3.65, 1.75
    rm = (r1 + r2) / 2
    cx = [0.65, 2.05, 3.45]
    ax.text(0.1, 6.6, "Radon-222 chain", fontsize=19, color=NAVY, fontweight="bold", va="center")
    node(cx[0], r0, 222, "Rn", "3.82 d", "alpha")
    node(cx[0], r1, 218, "Po", "3.10 min", "alpha")
    node(cx[0], r2, 214, "Pb", "26.8 min", "bkg")
    node(cx[1], rm, 214, "Bi", "19.9 min", "beta")
    node(cx[2], r1, 214, "Po", r"164 $\mu$s", "alpha")
    node(cx[2], r2, 210, "Pb", "22.2 yr", "end")
    adown(cx[0], r0, r1, r"$\alpha$ 5.59 MeV")
    adown(cx[0], r1, r2, r"$\alpha$ 6.11 MeV")
    bup(cx[0], r2, cx[1], rm)
    bup(cx[1], rm, cx[2], r1)
    adown(cx[2], r1, r2, r"$\alpha$ 7.83 MeV")
    ax.text(cx[0] - BW / 2, r2 - 0.5, "low-energy\nbackground", ha="left", va="top", fontsize=FS_L,
            color=ORANGE, fontweight="bold", linespacing=1.05)
    ax.text((cx[1] + cx[2]) / 2 + 0.15, r0, "BiPo: " + r"$^{214}$Po" + "\nfollows within\n164 " + r"$\mu$s",
            ha="center", va="center", fontsize=FS_T, color=NAVY, linespacing=1.15)

    dx = [5.35, 6.75, 8.15]
    ax.text(4.8, 6.6, "Radon-220 chain", fontsize=19, color=NAVY, fontweight="bold", va="center")
    node(dx[0], r0, 220, "Rn", "55.6 s", "alpha")
    node(dx[0], r1, 216, "Po", "0.145 s", "alpha")
    node(dx[0], r2, 212, "Pb", "10.6 h", "beta")
    node(dx[1], rm, 212, "Bi", "60.6 min", "beta")
    node(dx[2], r1, 212, "Po", r"0.3 $\mu$s", "beta")
    node(dx[1], 0.95, 208, "Tl", "3.05 min", "end")
    adown(dx[0], r0, r1, r"$\alpha$ 6.40 MeV")
    adown(dx[0], r1, r2, r"$\alpha$ 6.91 MeV")
    bup(dx[0], r2, dx[1], rm)
    bup(dx[1], rm, dx[2], r1, "β 64%")
    arrow(ax, dx[1], rm - BH / 2 - 0.04, dx[1], 0.95 + BH / 2 + 0.05, color=GREY, lw=1.8, ms=14)
    ax.text(dx[1] + 0.1, (rm + 0.95) / 2, r"$\alpha$ 36%", ha="left", va="center", fontsize=FS_T, color=GREY)
    ax.text(dx[2], r1 - 0.5, "8.95 MeV " + r"$\alpha$" + ",\nmerges with " + r"$\beta$", ha="center", va="top",
            fontsize=FS_T - 1, color=GREY, linespacing=1.1)
    ax.text((dx[1] + dx[2]) / 2 + 0.15, r0, "stable\nreference", ha="center", va="center",
            fontsize=FS_T, color=NAVY, linespacing=1.15)
    fig.savefig(OUT + "decay_chain_poster.png", dpi=300, facecolor="white")
    plt.close(fig)


def tpc_poster(width=9.2):
    H = 6.6
    fig = plt.figure(figsize=(width, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    FS = 17
    x0, x1 = 1.15, 4.75
    yb, yt = 0.2, 6.4
    y_cath, y_gate, y_liq, y_anode = 0.75, 5.0, 5.22, 5.45
    ax.add_patch(Rectangle((x0, yb), x1 - x0, yt - yb, fc="white", ec=NAVY, lw=2.2))
    ax.add_patch(Rectangle((x0, y_cath - 0.05), x1 - x0, y_liq - y_cath + 0.05, fc="#DCE9F7", ec="none"))
    ax.add_patch(Rectangle((x0, y_liq), x1 - x0, yt - 0.48 - y_liq, fc="#F4F8FC", ec="none"))
    ax.plot([x0, x1], [y_liq, y_liq], color="#3D78B8", lw=2)
    for yy in (yt - 0.26, yb + 0.26):
        for xc in np.linspace(x0 + 0.3, x1 - 0.3, 7):
            ax.add_patch(Circle((xc, yy), 0.2, fc=GOLD, ec=NAVY, lw=1.1))
    for yy, lab in ((y_anode, "anode"), (y_gate, "gate"), (y_cath, "cathode")):
        ax.plot([x0, x1], [yy, yy], color=GREY, lw=1.5, ls=(0, (2, 2)))
        ax.text(x0 - 0.1, yy + (0.1 if lab == "anode" else (-0.1 if lab == "gate" else 0)), lab,
                ha="right", va="center", fontsize=FS - 2, color=GREY)
    ax.text(x0 - 0.1, (y_cath + y_gate) / 2 + 0.6, "liquid\nxenon", ha="right", va="center", fontsize=FS - 2, color=NAVY)
    ax.text(x0 - 0.1, y_anode + 0.42, "gas", ha="right", va="center", fontsize=FS - 2, color=NAVY)
    # field labels
    ax.text(x1 + 0.12, (y_gate + y_anode) / 2, "extraction\nfield", ha="left", va="center", fontsize=FS - 3,
            color=ORANGE, linespacing=1.0)
    ax.plot([x1 + 0.05, x1 + 0.08, x1 + 0.08, x1 + 0.05], [y_gate, y_gate, y_anode, y_anode], color=ORANGE, lw=1.4)
    arrow(ax, x0 + 0.25, y_gate - 0.2, x0 + 0.25, y_cath + 0.25, color=GREY, lw=1.4, ms=12)
    ax.text(x0 + 0.36, 2.2, "drift field", rotation=90, fontsize=FS - 4, color=GREY, va="center")
    # interaction
    ix, iy = 3.0, 1.9
    star = Polygon([(ix + 0.27 * np.cos(tt) * (1 if k % 2 == 0 else 0.45),
                     iy + 0.27 * np.sin(tt) * (1 if k % 2 == 0 else 0.45))
                    for k, tt in enumerate(np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, 11)[:-1])],
                   fc=ORANGE, ec=NAVY, lw=1)
    ax.add_patch(star)
    for tx in (1.5, 2.3, 4.25):
        arrow(ax, ix, iy + 0.27, tx, yt - 0.48, color=TEAL, lw=1.3, ms=11)
    for bx in (2.0, 3.5, 4.3):
        arrow(ax, ix, iy - 0.27, bx, yb + 0.48, color=TEAL, lw=1.3, ms=11)
    ax.text(ix + 0.32, iy - 0.08, "S1", fontsize=FS + 3, color=TEAL, fontweight="bold", va="center")
    arrow(ax, ix, iy + 0.3, ix, y_gate - 0.04, color="#3D78B8", lw=1.8, ms=13, ls=(0, (4, 3)))
    ax.text(ix + 0.1, 3.5, r"e$^-$ drift", rotation=90, fontsize=FS - 3, color="#3D78B8", ha="left", va="center")
    ax.add_patch(Circle((ix, (y_liq + y_anode) / 2 + 0.03), 0.17, fc="#FFE9A8", ec=ORANGE, lw=1.4))
    ax.text(ix + 0.25, y_anode + 0.25, "S2", fontsize=FS + 3, color=ORANGE, fontweight="bold", va="center")

    # right: waveform + AFT
    wx0, wx1, wy0 = 5.75, width - 0.1, 4.1
    ax.text(wx0, 6.15, "Waveform", fontsize=FS + 1, color=NAVY, fontweight="bold", va="center")
    tt = np.linspace(0, 1, 600)
    sig = 0.45 * np.exp(-0.5 * ((tt - 0.1) / 0.01) ** 2) + 1.35 * np.exp(-0.5 * ((tt - 0.78) / 0.045) ** 2)
    ax.plot(wx0 + tt * (wx1 - wx0), wy0 + sig, color=NAVY, lw=1.6)
    ax.plot([wx0, wx1], [wy0, wy0], color=GREY, lw=1)
    ax.text(wx0 + 0.1 * (wx1 - wx0), wy0 + 0.6, "S1", ha="center", fontsize=FS, color=TEAL, fontweight="bold")
    ax.text(wx0 + 0.78 * (wx1 - wx0) - 0.3, wy0 + 1.25, "S2", ha="right", fontsize=FS, color=ORANGE, fontweight="bold")
    arrow(ax, wx0 + 0.1 * (wx1 - wx0), wy0 - 0.2, wx0 + 0.78 * (wx1 - wx0), wy0 - 0.2,
          color=GREY, lw=1.2, style="<|-|>", ms=9)
    ax.text(wx0 + 0.44 * (wx1 - wx0), wy0 - 0.45, "drift time  →  depth", ha="center", va="top", fontsize=FS - 3, color=GREY)

    box(ax, (wx0 + wx1) / 2, 1.65, wx1 - wx0 + 0.1, 2.5, GREY_L, NAVY, lw=1.3, r=0.08)
    ax.text(wx0 + 0.1, 2.55, "AFT", fontsize=FS + 1, color=NAVY, fontweight="bold", va="center")
    ax.text(wx0 + 0.1, 1.95, "= S1 on top array\n   / total S1", fontsize=FS - 2, color=NAVY, va="center", linespacing=1.15)
    ax.text(wx0 + 0.1, 0.95, "higher near the top:\ndepth from S1 alone", fontsize=FS - 3, color=GREY,
            va="center", linespacing=1.15)
    fig.savefig(OUT + "tpc_schematic_poster.png", dpi=300, facecolor="white")
    plt.close(fig)


# ======================================================================
# v2 (two-column A0 poster): horizontal decay chain + compact TPC
# ======================================================================
def decay_chain_h(width=14.4):
    Hh = 5.3
    fig = plt.figure(figsize=(width, Hh))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width); ax.set_ylim(0, Hh); ax.set_aspect("equal"); ax.axis("off")
    BW, BH = 1.3, 0.9
    FI, FT, FL = 26, 17, 19
    xs = [0.7, 3.3, 5.9, 8.5, 11.1, 13.7]

    def node(x, y, a, el, tt, kind):
        fc, ec = {"alpha": (TEAL_L, TEAL), "bkg": (ORANGE_L, ORANGE),
                  "beta": (GREY_L, GREY), "end": ("white", GREY)}[kind]
        box(ax, x, y, BW, BH, fc, ec, lw=2.6 if kind in ("alpha", "bkg") else 1.6, r=0.07)
        ax.text(x, y + 0.14, iso(a, el), ha="center", va="center", fontsize=FI, color=NAVY, fontweight="bold")
        ax.text(x, y - 0.25, tt, ha="center", va="center", fontsize=FT - 1, color=GREY)

    def step(x0, x1, y, lab, a=True):
        c = TEAL if a else GREY
        arrow(ax, x0 + BW / 2 + 0.05, y, x1 - BW / 2 - 0.07, y, color=c, lw=2.4, ms=18)
        ax.text((x0 + x1) / 2, y + 0.12, lab, ha="center", va="bottom", fontsize=FL, color=c,
                fontweight="bold" if a else "normal", linespacing=1.0)

    y1, y2 = 3.95, 1.45
    ax.text(0.08, 5.1, "Radon-222 chain", fontsize=22, color=NAVY, fontweight="bold", va="center")
    node(xs[0], y1, 222, "Rn", "3.82 d", "alpha")
    node(xs[1], y1, 218, "Po", "3.10 min", "alpha")
    node(xs[2], y1, 214, "Pb", "26.8 min", "bkg")
    node(xs[3], y1, 214, "Bi", "19.9 min", "beta")
    node(xs[4], y1, 214, "Po", r"164 $\mu$s", "alpha")
    node(xs[5], y1, 210, "Pb", "22.2 yr", "end")
    step(xs[0], xs[1], y1, r"$\alpha$" + "\n5.59 MeV")
    step(xs[1], xs[2], y1, r"$\alpha$" + "\n6.11 MeV")
    step(xs[2], xs[3], y1, r"$\beta$", a=False)
    step(xs[3], xs[4], y1, r"$\beta$", a=False)
    step(xs[4], xs[5], y1, r"$\alpha$" + "\n7.83 MeV")
    ax.text(xs[2], y1 - BH / 2 - 0.12, "low-energy background", ha="center", va="top", fontsize=FL,
            color=ORANGE, fontweight="bold")
    bx0, bx1, by = xs[3] - BW / 2, xs[4] + BW / 2, y1 - BH / 2 - 0.15
    ax.plot([bx0, bx0, bx1, bx1], [by + 0.08, by, by, by + 0.08], color=NAVY, lw=1.5)
    ax.text((bx0 + bx1) / 2, by - 0.08, "BiPo: two decays within ~0.2 ms", ha="center", va="top",
            fontsize=FT, color=NAVY)

    ax.text(0.08, 2.55, "Radon-220 chain (thoron)", fontsize=22, color=NAVY, fontweight="bold", va="center")
    node(xs[0], y2, 220, "Rn", "55.6 s", "alpha")
    node(xs[1], y2, 216, "Po", "0.145 s", "alpha")
    node(xs[2], y2, 212, "Pb", "10.6 h", "beta")
    node(xs[3], y2, 212, "Bi", "60.6 min", "beta")
    step(xs[0], xs[1], y2, r"$\alpha$" + "\n6.40 MeV")
    step(xs[1], xs[2], y2, r"$\alpha$" + "\n6.91 MeV")
    step(xs[2], xs[3], y2, r"$\beta$", a=False)
    bxs = xs[3] + BW / 2 + 0.05
    arrow(ax, bxs, y2 + 0.15, bxs + 0.75, y2 + 0.55, color=GREY, lw=1.8, ms=14)
    arrow(ax, bxs, y2 - 0.15, bxs + 0.75, y2 - 0.55, color=GREY, lw=1.8, ms=14)
    ax.text(bxs + 0.85, y2 + 0.6, r"$\alpha$ 36%  $\rightarrow$ $^{208}$Tl", va="center", fontsize=FT, color=GREY)
    ax.text(bxs + 0.85, y2 - 0.6, r"$\beta$ 64%  $\rightarrow$ $^{212}$Po ($\alpha$ 8.95 MeV, 0.3 $\mu$s)",
            va="center", fontsize=FT, color=GREY)

    ly = 0.25
    for x, fc, ec, lab in [(0.1, TEAL_L, TEAL, "alpha peak used in the fit"),
                           (5.0, ORANGE_L, ORANGE, "background of interest"),
                           (9.6, GREY_L, GREY, "not fitted")]:
        box(ax, x + 0.25, ly, 0.5, 0.34, fc, ec, lw=1.8, r=0.05)
        ax.text(x + 0.65, ly, lab, va="center", fontsize=FT, color=NAVY)
    fig.savefig(OUT + "decay_chain_h.png", dpi=300, facecolor="white")
    plt.close(fig)


def tpc_compact(width=6.9):
    H = 6.0
    fig = plt.figure(figsize=(width, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width); ax.set_ylim(0, H); ax.set_aspect("equal"); ax.axis("off")
    FS = 17
    x0, x1 = 0.95, 3.85
    yb, yt = 0.1, 5.9
    y_cath, y_gate, y_liq, y_anode = 0.62, 4.55, 4.76, 4.97
    ax.add_patch(Rectangle((x0, yb), x1 - x0, yt - yb, fc="white", ec=NAVY, lw=2.2))
    ax.add_patch(Rectangle((x0, y_cath - 0.05), x1 - x0, y_liq - y_cath + 0.05, fc="#DCE9F7", ec="none"))
    ax.add_patch(Rectangle((x0, y_liq), x1 - x0, yt - 0.44 - y_liq, fc="#F4F8FC", ec="none"))
    ax.plot([x0, x1], [y_liq, y_liq], color="#3D78B8", lw=2)
    for yy in (yt - 0.24, yb + 0.24):
        for xc in np.linspace(x0 + 0.27, x1 - 0.27, 6):
            ax.add_patch(Circle((xc, yy), 0.19, fc=GOLD, ec=NAVY, lw=1.1))
    for yy, lab, dy in ((y_anode, "anode", 0.12), (y_gate, "gate", -0.12), (y_cath, "cathode", 0)):
        ax.plot([x0, x1], [yy, yy], color=GREY, lw=1.4, ls=(0, (2, 2)))
        ax.text(x0 - 0.08, yy + dy, lab, ha="right", va="center", fontsize=FS - 3, color=GREY)
    ax.text(x0 - 0.08, 2.9, "liquid\nxenon", ha="right", va="center", fontsize=FS - 2, color=NAVY)
    ax.text(x0 - 0.08, y_anode + 0.42, "gas", ha="right", va="center", fontsize=FS - 2, color=NAVY)
    # electric field direction (points down; electrons move up)
    arrow(ax, x1 + 0.13, y_anode + 0.02, x1 + 0.13, y_gate - 0.02, color=ORANGE, lw=1.8, ms=11)
    ax.text(x1 + 0.25, (y_gate + y_anode) / 2, "extraction\nfield", ha="left", va="center",
            fontsize=FS - 4, color=ORANGE, linespacing=1.0)
    arrow(ax, x1 + 0.13, y_gate - 0.12, x1 + 0.13, y_cath + 0.02, color=GREY, lw=1.8, ms=12)
    ax.text(x1 + 0.25, 2.6, "drift\nfield", ha="left", va="center", fontsize=FS - 4, color=GREY, linespacing=1.0)
    ix, iy = 2.5, 1.7
    # explosion icon
    radii_o = [0.34, 0.24, 0.38, 0.27, 0.33, 0.22, 0.36, 0.25, 0.31, 0.23, 0.37, 0.26]
    angs = np.linspace(0, 2 * np.pi, 25)[:-1] + 0.1
    pts = [(ix + (radii_o[k // 2] if k % 2 == 0 else 0.13) * np.cos(a),
            iy + (radii_o[k // 2] if k % 2 == 0 else 0.13) * np.sin(a)) for k, a in enumerate(angs)]
    ax.add_patch(Polygon(pts, fc=ORANGE, ec="#8A3A0F", lw=1, zorder=5))
    pts2 = [(ix + (0.6 * radii_o[k // 2] if k % 2 == 0 else 0.07) * np.cos(a + 0.13),
             iy + (0.6 * radii_o[k // 2] if k % 2 == 0 else 0.07) * np.sin(a + 0.13)) for k, a in enumerate(angs)]
    ax.add_patch(Polygon(pts2, fc="#FFD23F", ec="none", zorder=6))
    for tx in (1.25, 1.95, 3.55):
        arrow(ax, ix, iy + 0.25, tx, yt - 0.44, color=TEAL, lw=1.2, ms=10)
    for bx in (1.6, 3.0, 3.6):
        arrow(ax, ix, iy - 0.25, bx, yb + 0.44, color=TEAL, lw=1.2, ms=10)
    ax.text(ix + 0.4, iy - 0.05, "S1", fontsize=FS + 3, color=TEAL, fontweight="bold", va="center")
    arrow(ax, ix, iy + 0.36, ix, y_liq - 0.2, color="#3D78B8", lw=1.8, ms=12, ls=(0, (4, 3)))
    ax.text(ix + 0.08, 3.2, r"e$^-$ drift", rotation=90, fontsize=FS - 4, color="#3D78B8", ha="left", va="center")
    ey = (y_liq + y_anode) / 2 + 0.02
    ax.add_patch(Circle((ix, ey), 0.1, fc="#3D78B8", ec=NAVY, lw=0.9, zorder=6))
    ax.text(ix, ey - 0.005, r"e$^-$", ha="center", va="center", fontsize=7.5, color="white", fontweight="bold", zorder=7)
    for a in np.linspace(0.35, np.pi - 0.35, 5):
        ax.plot([ix + 0.15 * np.cos(a), ix + 0.3 * np.cos(a)], [ey + 0.15 * np.sin(a), ey + 0.3 * np.sin(a)],
                color=ORANGE, lw=1.5, zorder=6)
    ax.text(ix + 0.25, y_anode + 0.28, "S2", fontsize=FS + 3, color=ORANGE, fontweight="bold", va="center")
    # waveform (right, top)
    wx0, wx1, wy0 = 5.0, width - 0.05, 4.2
    tt = np.linspace(0, 1, 500)
    sig = 0.4 * np.exp(-0.5 * ((tt - 0.1) / 0.012) ** 2) + 1.15 * np.exp(-0.5 * ((tt - 0.78) / 0.05) ** 2)
    ax.plot(wx0 + tt * (wx1 - wx0), wy0 + sig, color=NAVY, lw=1.6)
    ax.plot([wx0, wx1], [wy0, wy0], color=GREY, lw=1)
    ax.text(wx0 + 0.1 * (wx1 - wx0), wy0 + 0.5, "S1", ha="center", fontsize=FS - 1, color=TEAL, fontweight="bold")
    ax.text(wx0 + 0.78 * (wx1 - wx0) - 0.22, wy0 + 1.0, "S2", ha="right", fontsize=FS - 1, color=ORANGE, fontweight="bold")
    arrow(ax, wx0 + 0.1 * (wx1 - wx0), wy0 - 0.18, wx0 + 0.78 * (wx1 - wx0), wy0 - 0.18, color=GREY, lw=1.1, style="<|-|>", ms=8)
    ax.text((wx0 + wx1) / 2, wy0 - 0.38, "drift time " + r"$\propto$" + " depth", ha="center", va="top", fontsize=FS - 4, color=GREY)
    # compact AFT box (right, bottom), sized to its text
    bx0, bw, by0, bh = wx0 - 0.1, wx1 - wx0 + 0.15, 1.85, 1.15
    ax.add_patch(FancyBboxPatch((bx0, by0), bw, bh, boxstyle="round,pad=0,rounding_size=0.08",
                                fc=GREY_L, ec=NAVY, lw=1.2))
    ax.text(bx0 + bw / 2, by0 + bh / 2, r"AFT $=\dfrac{\mathrm{PE_{top}}}{\mathrm{PE_{total}}}$", ha="center", va="center",
            fontsize=FS, color=NAVY)
    ax.text(bx0 + bw / 2, by0 - 0.12, "higher AFT " + r"$\rightarrow$" + " closer\nto the top: indicates depth",
            ha="center", va="top", fontsize=FS - 5, color=GREY, linespacing=1.1)
    ax.text(bx0 + bw / 2, by0 - 0.75, "PE: photoelectrons", ha="center", va="top", fontsize=FS - 6, color=GREY)
    fig.savefig(OUT + "tpc_compact.png", dpi=300, facecolor="white")
    plt.close(fig)
