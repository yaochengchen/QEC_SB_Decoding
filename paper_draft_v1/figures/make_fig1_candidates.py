"""Three schematic, publication-oriented layouts for manuscript Figure 1.

All depicted masses are illustrative. No experimental values are encoded.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


OUT = Path(__file__).with_name("figure1_candidates")
OUT.mkdir(exist_ok=True)

INK = "#203145"
MUTED = "#597086"
BLUE = "#2874B8"
BLUE_LIGHT = "#EAF3FB"
TEAL = "#188C83"
TEAL_LIGHT = "#E9F7F3"
GOLD = "#D29329"
GOLD_LIGHT = "#FFF5E3"
PURPLE = "#7955A5"
PURPLE_LIGHT = "#F3EEFA"
RED = "#C95559"
RED_LIGHT = "#FCEEEE"
PALE = "#F5F8FB"
LINE = "#D8E2EC"
WHITE = "#FFFFFF"
COLORS = [BLUE, TEAL, GOLD, PURPLE]
PALES = [BLUE_LIGHT, TEAL_LIGHT, GOLD_LIGHT, PURPLE_LIGHT]
MASS = [0.48, 0.83, 0.59, 0.68]


def canvas():
    fig, ax = plt.subplots(figsize=(16, 8.5), dpi=180)
    fig.patch.set_facecolor(WHITE)
    ax.set(xlim=(0, 100), ylim=(0, 58))
    ax.set_aspect("equal")
    ax.axis("off")
    fig.subplots_adjust(left=0.015, right=0.985, bottom=0.03, top=0.98)
    return fig, ax


def box(ax, x, y, w, h, fc=WHITE, ec=LINE, lw=1.2, radius=1.3, z=1):
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
    )
    ax.add_patch(p)
    return p


def txt(ax, x, y, s, size=11, color=INK, weight="normal", ha="left", va="center", **kw):
    return ax.text(x, y, s, fontsize=size, color=color, weight=weight,
                   ha=ha, va=va, family="DejaVu Sans", zorder=10, **kw)


def arrow(ax, x1, y1, x2, y2, color=MUTED, lw=1.6, scale=13, z=5):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=scale, linewidth=lw, color=color,
                                 shrinkA=0, shrinkB=0, zorder=z))


def dot(ax, x, y, r=0.62, fc=WHITE, ec=BLUE, lw=1.2):
    ax.add_patch(Circle((x, y), r, facecolor=fc, edgecolor=ec, linewidth=lw, zorder=6))


def patch_icon(ax, cx, cy, scale=1, c=BLUE):
    """A generic rotated-code motif; red marks denote a measured syndrome."""
    pts = []
    for i in range(4):
        for j in range(4):
            x = cx + (i - 1.5) * 2.35 * scale
            y = cy + (j - 1.5) * 2.35 * scale
            pts.append((x, y))
            if i < 3:
                ax.plot([x, x + 2.35 * scale], [y, y], color="#ACC5DC", lw=1.0, zorder=2)
            if j < 3:
                ax.plot([x, x], [y, y + 2.35 * scale], color="#ACC5DC", lw=1.0, zorder=2)
    for x, y in pts:
        dot(ax, x, y, r=0.37 * scale, fc=WHITE, ec=c, lw=1.0)
    for dx, dy in [(-2.35, 0), (2.35, 2.35)]:
        ax.add_patch(Circle((cx + dx * scale, cy + dy * scale), 0.73 * scale,
                            facecolor=RED, edgecolor=WHITE, linewidth=0.8, zorder=7))


def candidate_icon(ax, x, y, c, compact=False):
    offsets = [(-2.4, 0.8), (-0.6, -0.65), (1.1, 0.95), (2.8, -0.25)] if not compact else [(-1.6, .4), (0, -.5), (1.6, .5)]
    for k, (dx, dy) in enumerate(offsets):
        dot(ax, x + dx, y + dy, r=0.5 if compact else 0.6,
            fc=c if k == 1 else WHITE, ec=c, lw=1.3)


def orbit_icon(ax, x, y, c, compact=False):
    s = 0.74 if compact else 1.0
    for dx, dy in [(-3, 0), (-1.8, 1.9), (0, 0), (1.8, 1.9), (3, 0), (-1.8, -1.9), (1.8, -1.9)]:
        dot(ax, x + dx * s, y + dy * s, r=0.49 * s,
            fc=c if dx == 0 and dy == 0 else WHITE, ec=c, lw=1.15)
    ax.add_patch(Circle((x, y), 3.6 * s, facecolor="none", edgecolor=c,
                        linewidth=1.05, linestyle=(0, (3, 3)), alpha=.75, zorder=4))


def header(ax, letter, title, subtitle):
    heading = f"{letter}  {title}" if letter else title
    txt(ax, 3.2, 54.7, heading, 20, INK, "bold")
    txt(ax, 3.2, 51.5, subtitle, 10.5, MUTED)
    ax.plot([3.2, 96.8], [49.6, 49.6], color=LINE, lw=1.1)


def footer(ax):
    ax.plot([3.2, 96.8], [5.1, 5.1], color=LINE, lw=1.0)
    txt(ax, 3.2, 2.5,
        "Schematic only: bar lengths are illustrative; the radius is measured in generator coordinates.",
        9.2, MUTED)


def save(fig, stem, dest=OUT):
    fig.savefig(dest / f"{stem}.png", dpi=220, facecolor=WHITE)
    fig.savefig(dest / f"{stem}.svg", facecolor=WHITE)
    plt.close(fig)


def variant_1():
    """Editorial pipeline: patch -> four lanes -> mass comparison."""
    fig, ax = canvas()
    header(ax, "01", "Four parallel logical-class searches",
           "A single measured syndrome defines four feasible sectors; each sector is searched and locally scored.")

    box(ax, 3, 10, 16, 36, PALE, LINE, 1.3)
    txt(ax, 11, 42.9, "Measured syndrome", 12, INK, "bold", "center")
    patch_icon(ax, 11, 29.7, 1.0)
    txt(ax, 11, 18.0, "$s$  +  feasible seed $e_0$", 10.6, INK, ha="center")
    txt(ax, 11, 14.6, "code-capacity noise", 9.3, MUTED, ha="center")
    arrow(ax, 19.0, 28, 23.1, 28, BLUE, 2.0)

    txt(ax, 25.0, 45.5, "FEASIBLE CLASS", 9.5, MUTED, "bold")
    txt(ax, 40.0, 45.5, "SB TERMINAL POOL", 9.5, MUTED, "bold")
    txt(ax, 60.8, 45.5, "LOCAL ORBIT  $d_H(u,u_L^*) \\leq 2$", 9.5, MUTED, "bold")
    txt(ax, 81.5, 45.5, "LOCAL MASS", 9.5, MUTED, "bold")
    for i, yy in enumerate([39.7, 32.4, 25.1, 17.8]):
        c, pale = COLORS[i], PALES[i]
        box(ax, 23.5, yy-3.0, 73.2, 6.05, pale, c, .95, .75)
        box(ax, 24.5, yy-2.05, 11, 4.1, WHITE, c, 1.1, .7)
        txt(ax, 30, yy, f"$L={i}$", 13.3, c, "bold", "center")
        arrow(ax, 36.2, yy, 39, yy, c, 1.25, 10)
        candidate_icon(ax, 45.0, yy, c, compact=True)
        txt(ax, 52.3, yy, "$u_L^*$", 11.5, c, "bold")
        arrow(ax, 55.7, yy, 59.4, yy, c, 1.25, 10)
        orbit_icon(ax, 65.0, yy, c, compact=True)
        arrow(ax, 70.1, yy, 75, yy, c, 1.25, 10)
        box(ax, 76.3, yy-1.25, 17.8, 2.5, WHITE, LINE, .6, .35)
        ax.add_patch(Rectangle((76.6, yy-1.0), MASS[i]*15.7, 2.0,
                               facecolor=c, edgecolor="none", zorder=3))
        if i == 1:
            txt(ax, 95.2, yy, "max", 9.3, c, "bold")

    box(ax, 23.5, 10, 73.2, 3.85, WHITE, LINE, 1.0, .65)
    txt(ax, 25.2, 11.92, "Binary states preserve $s$ and $L$; return the winning center as correction.", 10.0, INK)
    txt(ax, 78.0, 11.92, "$\\widehat L=\\arg\\max_L \\widehat Z_L^{(2)}$", 10.5, BLUE, "bold")
    footer(ax)
    save(fig, "fig1_A1_parallel_lanes")


def variant_2():
    """Transit map: a compact process route repeated in four color-coded lanes."""
    fig, ax = canvas()
    header(ax, "02", "From one syndrome to one correction",
           "A repeated route makes the four class-wise computations and the final decision visible at a glance.")
    box(ax, 3.2, 38.2, 93.6, 8.3, BLUE_LIGHT, BLUE, 1.3, 1.1)
    patch_icon(ax, 9.4, 42.4, .54)
    txt(ax, 17, 42.6, "Measure $s$", 13, INK, "bold")
    arrow(ax, 28.8, 42.4, 32.5, 42.4, BLUE, 1.9)
    txt(ax, 52, 42.4, "Construct four syndrome-compatible class seeds", 11.6, INK, ha="center")
    arrow(ax, 72.5, 42.4, 77, 42.4, BLUE, 1.9)
    box(ax, 77.5, 39.8, 15.9, 5.2, WHITE, BLUE, 1, .75)
    txt(ax, 85.5, 42.4, "$e_L(u)=e_L(0)\\oplus uG$", 10.5, BLUE, ha="center")

    xcols = [12.3, 31.2, 51.0, 69.5, 87.8]
    labels = ["CLASS", "SB TERMINAL STATES", "LOWEST-WEIGHT CENTER", "RADIUS-2 SOE", "LOCAL MASS"]
    for x, label in zip(xcols, labels):
        txt(ax, x, 35.1, label, 9.2, MUTED, "bold", "center")
    for i, yy in enumerate([31.0, 25.0, 19.0, 13.0]):
        c = COLORS[i]
        ax.plot([18.7, 82.2], [yy, yy], color=c, lw=2.25, alpha=.67, zorder=2)
        for x in [21.5, 40.6, 60.2, 79.8]:
            dot(ax, x, yy, .30, WHITE, c, 1.2)
        box(ax, 5.5, yy-2.0, 13.2, 4, PALES[i], c, 1.2, 1.1)
        txt(ax, 12.1, yy, f"$L={i}$", 13.2, c, "bold", "center")
        candidate_icon(ax, 31.2, yy, c, compact=True)
        dot(ax, 51.0, yy, .95, c, WHITE, 1.5)
        txt(ax, 51.0, yy-.05, "*", 12, WHITE, "bold", "center")
        orbit_icon(ax, 69.5, yy, c, compact=True)
        box(ax, 82.3, yy-1.25, 12.4, 2.5, WHITE, LINE, .8, .3)
        ax.add_patch(Rectangle((82.6, yy-1.0), MASS[i]*10.5, 2.0, facecolor=c, edgecolor="none", zorder=3))
        if i == 1:
            ax.add_patch(Circle((95.5, yy), 1.03, fc=TEAL, ec=WHITE, lw=1, zorder=6))
            txt(ax, 95.5, yy, "✓", 10.4, WHITE, "bold", "center")
    box(ax, 5.5, 6.5, 89.2, 3.6, PALE, LINE, 1.0, .7)
    txt(ax, 8, 8.25, "Read each lane:", 10.2, INK, "bold")
    txt(ax, 21, 8.25, "search terminal pool", 10.2, MUTED)
    arrow(ax, 39, 8.25, 42.3, 8.25, MUTED, 1.1, 9)
    txt(ax, 44, 8.25, "select $u_L^*$", 10.2, MUTED)
    arrow(ax, 57, 8.25, 60.3, 8.25, MUTED, 1.1, 9)
    txt(ax, 62, 8.25, "sum $w^{|e_L(u)|}$", 10.2, MUTED)
    arrow(ax, 76, 8.25, 79.3, 8.25, MUTED, 1.1, 9)
    txt(ax, 80.4, 8.25, "return winning center", 9.6, MUTED)
    footer(ax)
    save(fig, "fig1_A2_transit_map")


def variant_3(final=False):
    """Magazine layout: large central class matrix plus a decisive output block."""
    fig, ax = canvas()
    header(ax, "" if final else "03", "Search first, then aggregate",
           "Each class yields one searched center; nearby stabilizer-equivalent errors contribute to its local score.")
    box(ax, 3.1, 10.1, 19.0, 36.7, PALE, LINE, 1.15)
    txt(ax, 12.6, 43.4, "INPUT", 10, BLUE, "bold", "center")
    patch_icon(ax, 12.6, 32.4, 1.03)
    txt(ax, 12.6, 20.7, "one syndrome $s$", 11.6, INK, "bold", "center")
    txt(ax, 12.6, 16.7, "feasible seed $e_0$", 10.2, MUTED, ha="center")
    arrow(ax, 22.2, 28.3, 26.3, 28.3, BLUE, 1.7)

    box(ax, 26.8, 10.1, 47.7, 36.7, WHITE, LINE, 1.2)
    txt(ax, 29.1, 43.4, "FOUR INDEPENDENT LOGICAL-CLASS LANES", 10, BLUE, "bold")
    txt(ax, 39, 39.8, "CP-SB candidate pool", 9.1, MUTED, "bold", "center")
    txt(ax, 55.2, 39.8, "selected center", 9.1, MUTED, "bold", "center")
    txt(ax, 67.6, 39.8, "local sum ($\\leq 2$ flips)", 8.7, MUTED, "bold", "center")
    for i, yy in enumerate([35.0, 28.3, 21.6, 14.9]):
        c = COLORS[i]
        box(ax, 28.5, yy-2.6, 44.3, 5.2, PALES[i], "none", 0, .6)
        box(ax, 29.3, yy-1.6, 7.0, 3.2, WHITE, c, 1.0, .6)
        txt(ax, 32.8, yy, f"$L={i}$", 11.8, c, "bold", "center")
        candidate_icon(ax, 43.5, yy, c, compact=True)
        arrow(ax, 48.1, yy, 51.4, yy, c, 1.1, 9)
        dot(ax, 55.1, yy, .78, c, WHITE, 1.2)
        arrow(ax, 59.0, yy, 62.3, yy, c, 1.1, 9)
        orbit_icon(ax, 67.0, yy, c, compact=True)
    arrow(ax, 74.5, 28.3, 78.3, 28.3, BLUE, 1.7)

    box(ax, 78.9, 10.1, 17.9, 36.7, BLUE_LIGHT, BLUE, 1.25)
    txt(ax, 87.9, 43.4, "DECIDE", 10, BLUE, "bold", "center")
    txt(ax, 87.9, 40.2, "local class mass", 10.2, INK, ha="center")
    for i, yy in enumerate([35.1, 30.1, 25.1, 20.1]):
        c = COLORS[i]
        txt(ax, 81.2, yy, f"$L={i}$", 9.9, c, "bold")
        box(ax, 85.4, yy-1.0, 9.4, 2.0, WHITE, LINE, .5, .3)
        ax.add_patch(Rectangle((85.7, yy-.75), MASS[i]*8.45, 1.5, facecolor=c, edgecolor="none", zorder=3))
    txt(ax, 87.8, 15.7, "$\\widehat L=\\arg\\max_L \\widehat Z_L^{(2)}$", 10.0, BLUE, "bold", "center")
    txt(ax, 87.8, 12.4, "return center as correction", 9.4, INK, ha="center")
    footer(ax)
    if final:
        save(fig, "fig1_decoder_pipeline", Path(__file__).parent)
    else:
        save(fig, "fig1_A3_editorial_matrix")


if __name__ == "__main__":
    variant_1()
    variant_2()
    variant_3()
    variant_3(final=True)
    print("Wrote three PNG previews and three SVG source files to", OUT)
