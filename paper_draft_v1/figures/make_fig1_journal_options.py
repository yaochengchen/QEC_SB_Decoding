"""Two journal-style alternatives for manuscript Figure 1.

The drawings are schematic method diagrams and do not encode experimental data.
Regular labels remain editable in SVG output; mathtext may be exported as paths.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle


HERE = Path(__file__).parent
OUT = HERE / "figure1_journal_options"
OUT.mkdir(exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 7.5,
        "mathtext.fontset": "dejavusans",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)

BLACK = "#20252B"
GRAY = "#626B73"
LIGHT = "#C8CFD5"
FAINT = "#E7EAED"
BLUE = "#246A9B"
WHITE = "#FFFFFF"


def figure(width=7.1, height=3.6):
    fig, ax = plt.subplots(figsize=(width, height), dpi=180)
    fig.patch.set_facecolor(WHITE)
    ax.set_xlim(0, 710)
    ax.set_ylim(0, height * 100)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    return fig, ax


def label(ax, x, y, s, size=7.5, color=BLACK, weight="normal", ha="left", va="center"):
    return ax.text(x, y, s, fontsize=size, color=color, weight=weight, ha=ha, va=va, zorder=10)


def line(ax, x1, y1, x2, y2, color=LIGHT, lw=0.65, ls="-"):
    ax.plot([x1, x2], [y1, y2], color=color, linewidth=lw, linestyle=ls, zorder=1)


def arrow(ax, x1, y1, x2, y2, color=BLACK, lw=0.75):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=7.5,
            linewidth=lw, color=color, shrinkA=0, shrinkB=0, zorder=4
        )
    )


def circle(ax, x, y, r=3.4, edge=BLACK, fill=WHITE, lw=0.75):
    ax.add_patch(Circle((x, y), r, edgecolor=edge, facecolor=fill, linewidth=lw, zorder=5))


def patch(ax, cx, cy, size=11.5):
    """Minimal check lattice; blue disks denote schematic detection events."""
    for i in range(4):
        for j in range(4):
            x = cx + (i - 1.5) * size
            y = cy + (j - 1.5) * size
            if i < 3:
                line(ax, x, y, x + size, y, LIGHT, 0.55)
            if j < 3:
                line(ax, x, y, x, y + size, LIGHT, 0.55)
            circle(ax, x, y, 1.7, GRAY, WHITE, 0.55)
    for dx, dy in [(-size, 0), (size, size)]:
        circle(ax, cx + dx, cy + dy, 4.2, BLUE, BLUE, 0.6)


def candidate_pool(ax, cx, cy):
    for dx, dy in [(-13, 5), (-4, -4), (6, 4), (15, -3)]:
        circle(ax, cx + dx, cy + dy, 3.25, BLACK, WHITE, 0.7)


def local_set(ax, cx, cy):
    # One center and representative members at generator Hamming distances 1 and 2.
    circle(ax, cx - 12, cy, 3.5, BLUE, BLUE, 0.75)
    circle(ax, cx, cy + 5, 3.5, BLUE, WHITE, 0.75)
    circle(ax, cx + 12, cy, 3.5, BLUE, WHITE, 0.75)
    line(ax, cx - 8, cy + 1, cx - 4, cy + 4, BLUE, 0.7)
    line(ax, cx + 4, cy + 4, cx + 8, cy + 1, BLUE, 0.7)


def save(fig, name):
    fig.savefig(OUT / (name + ".png"), dpi=320, facecolor=WHITE)
    fig.savefig(OUT / (name + ".svg"), facecolor=WHITE)
    plt.close(fig)


def option_1():
    """Four equal lanes, as a compact journal line drawing."""
    fig, ax = figure(height=3.0)
    ax.set_ylim(45, 345)
    ytop = 308
    for x, s in [
        (19, "Input"),
        (167, "Logical class"),
        (276, "CP-SB candidates"),
        (375, "Center"),
        (460, "Local SOE"),
        (545, "Local mass"),
    ]:
        label(ax, x, ytop, s, 7.8, BLACK, "bold", "center")
    line(ax, 16, 294, 694, 294, BLACK, 0.8)

    patch(ax, 67, 177, 15.5)
    label(ax, 67, 259, "syndrome $s$", 8.1, BLACK, ha="center")
    label(ax, 67, 87, "feasible $e_0$", 8.0, GRAY, ha="center")
    arrow(ax, 107, 178, 133, 178, BLACK)

    rows = [259, 203, 147, 91]
    for idx, y in enumerate(rows):
        if idx:
            line(ax, 139, y + 28, 580, y + 28, FAINT, 0.55)
        label(ax, 166, y, f"$L={idx}$", 8.5, BLACK, ha="center")
        arrow(ax, 190, y, 217, y)
        candidate_pool(ax, 261, y)
        arrow(ax, 294, y, 322, y)
        circle(ax, 361, y, 4.3, BLACK, BLACK, 0.8)
        label(ax, 361, y - 15, "$u_L^*$", 7.5, BLACK, ha="center")
        arrow(ax, 386, y, 411, y)
        local_set(ax, 451, y)
        arrow(ax, 480, y, 504, y, BLUE)
        label(ax, 550, y, rf"$\widehat Z_{{{idx}}}^{{(2)}}$", 8.5, BLUE, ha="center")

    # Four symbolic scores feed a decision; no invented numerical bar lengths.
    line(ax, 580, 91, 580, 259, BLACK, 0.75)
    for y in rows:
        line(ax, 574, y, 580, y, BLACK, 0.75)
    arrow(ax, 580, 175, 601, 175)
    ax.add_patch(Rectangle((607, 120), 87, 109, fill=False, ec=BLACK, lw=0.8))
    label(ax, 650.5, 202, "select class", 8.2, BLACK, "bold", "center")
    label(ax, 650.5, 178, r"$\widehat L=\arg\max_L$", 7.4, BLACK, ha="center")
    label(ax, 650.5, 160, r"$\widehat Z_L^{(2)}$", 7.4, BLACK, ha="center")
    line(ax, 619, 145, 682, 145, LIGHT, 0.6)
    label(ax, 650.5, 132, "return correction", 7.3, GRAY, ha="center")
    save(fig, "fig1_J1_four_lanes")


def option_2():
    """Two compact panels: class construction and one repeated search rule."""
    fig, ax = figure(height=3.55)
    label(ax, 16, 328, "a", 9, BLACK, "bold")
    label(ax, 40, 328, "Four syndrome-compatible logical classes", 8.4, BLACK, "bold")
    line(ax, 16, 314, 694, 314, BLACK, 0.8)
    patch(ax, 73, 265, 10.5)
    label(ax, 73, 226, "measured $s$", 7.7, BLACK, ha="center")
    arrow(ax, 110, 265, 145, 265)
    label(ax, 159, 283, "feasible seed $e_0$", 7.8, GRAY)
    line(ax, 168, 265, 611, 265, LIGHT, 0.8)
    for idx, x in enumerate([224, 347, 470, 593]):
        circle(ax, x, 265, 4.0, BLACK, WHITE, 0.8)
        label(ax, x, 241, f"$L={idx}$", 8.4, BLACK, "bold", "center")
        label(ax, x, 222, rf"$e_{{{idx}}}(0)$", 7.6, GRAY, ha="center")
    label(ax, 682, 265, r"$s$ fixed", 7.4, BLUE, ha="right")

    label(ax, 16, 184, "b", 9, BLACK, "bold")
    label(ax, 40, 184, "Apply the same search and aggregation rule within each class", 8.4, BLACK, "bold")
    line(ax, 16, 170, 694, 170, BLACK, 0.8)
    label(ax, 72, 118, "CP-SB pool", 7.8, BLACK, "bold", "center")
    candidate_pool(ax, 72, 93)
    arrow(ax, 111, 93, 149, 93)
    label(ax, 199, 118, "select center", 7.8, BLACK, "bold", "center")
    circle(ax, 199, 93, 4.2, BLACK, BLACK, 0.8)
    label(ax, 199, 70, "$u_L^*$", 7.8, BLACK, ha="center")
    arrow(ax, 236, 93, 272, 93)
    label(ax, 337, 118, "radius-2 SOE", 7.8, BLACK, "bold", "center")
    local_set(ax, 337, 93)
    label(ax, 337, 69, "0, 1, or 2 generator flips", 7.2, GRAY, ha="center")
    arrow(ax, 392, 93, 429, 93, BLUE)
    label(ax, 483, 118, "local mass", 7.8, BLUE, "bold", "center")
    label(ax, 483, 93, r"$\widehat Z_L^{(2)}=\sum_{d_H(u,u_L^*)\leq2}w^{|e_L(u)|}$",
          7.1, BLUE, ha="center")
    arrow(ax, 545, 93, 574, 93)
    ax.add_patch(Rectangle((578, 63), 117, 60, fill=False, ec=BLACK, lw=0.8))
    label(ax, 636.5, 103, r"$\widehat L=\arg\max_L\widehat Z_L^{(2)}$",
          7.3, BLACK, ha="center")
    label(ax, 636.5, 80, "return center correction", 7.2, GRAY, ha="center")
    label(ax, 40, 25, r"Class-wise construction: $e_L(u)=e_L(0)\oplus uG$.",
          7.5, GRAY)
    save(fig, "fig1_J2_two_panels")


if __name__ == "__main__":
    option_1()
    option_2()
    print(OUT)
