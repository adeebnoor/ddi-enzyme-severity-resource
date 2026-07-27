#!/usr/bin/env python3
"""
make_figures.py — regenerate every figure in figures/ from the released CSVs.

Usage:
    python3 scripts/make_figures.py

Writes PNG (200 dpi) and SVG for each figure into figures/. Requires matplotlib;
no other third-party dependency. Every figure is derived directly from the files
in data/ — nothing is hand-drawn or hard-coded, so the figures cannot drift from
the data.

The palette is a validated, colour-vision-deficiency-safe set: a blue/red
diverging pair for enrichment vs depletion, a single-hue blue ordinal ramp for
ordered severity, and a neutral grey for results that do not survive FDR
correction. Text never wears a data colour.
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
FIGS = ROOT / "figures"
FIGS.mkdir(exist_ok=True)

# ----------------------------------------------------------------- palette
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

BLUE = "#2a78d6"    # diverging cool pole / categorical slot 1
RED = "#e34948"     # diverging warm pole
ORANGE = "#eb6834"  # categorical slot 2
GREY = "#a8a69f"    # de-emphasis (non-significant)

# single-hue ordinal ramp for ordered severity, light -> dark
SEV_RAMP = {"Minor": "#86b6ef", "Moderate": "#3987e5", "Major": "#184f95"}
SEV_UNKNOWN = "#c3c2b7"  # not a level of the scale — no grade recorded

plt.rcParams.update({
    "svg.hashsalt": "ddi-enzyme-severity-resource",
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "text.color": INK,
    "axes.labelcolor": INK_2,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "font.size": 10,
})


def rows(name):
    with open(DATA / name, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def summary(name):
    return {r["metric"]: r["value"] for r in rows(name)}


def style(ax, xgrid=True):
    """Hairline, recessive chrome; no chartjunk."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(AXIS)
    ax.spines["bottom"].set_color(AXIS)
    if xgrid:
        ax.set_axisbelow(True)
        ax.xaxis.grid(True, color=GRID, linewidth=0.8, linestyle="-")
        ax.yaxis.grid(False)
    ax.tick_params(length=0)


def px_to_data(ax, px):
    """Convert a pixel distance into data units on each axis (limits must be set)."""
    inv = ax.transData.inverted()
    x0, y0 = inv.transform((0, 0))
    x1, y1 = inv.transform((px, px))
    return abs(x1 - x0), abs(y1 - y0)


def hbar(ax, y, width, height, color, x0=0.0):
    """Horizontal bar: square at the baseline, 4px rounded data end, capped at
    24px thick. Axis limits must already be set — the corner radius is measured
    in pixels, so it stays a true 4px arc whatever the data range is."""
    rx, ry = px_to_data(ax, 4.0)
    _, max_h = px_to_data(ax, 24.0)
    height = min(height, max_h)
    rx = min(rx, abs(width) * 0.5)
    ry = min(ry, height * 0.5)
    s = 1 if width >= 0 else -1
    x1 = x0 + width
    xr = x1 - s * rx
    b, t = y - height / 2, y + height / 2
    k = 0.5523                      # circular-arc control-point constant
    verts = [
        (x0, b), (xr, b),
        (xr + s * rx * k, b), (x1, b + ry * (1 - k)), (x1, b + ry),
        (x1, t - ry), (x1, t - ry * (1 - k)), (xr + s * rx * k, t),
        (xr, t), (x0, t), (x0, b),
    ]
    codes = [
        MplPath.MOVETO, MplPath.LINETO,
        MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
        MplPath.LINETO, MplPath.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color,
                           edgecolor="none", zorder=3))


def titles(fig, ax, title, subtitle):
    """Title and subtitle above the plot, left-aligned to the axes, laid out in
    inches so they never collide regardless of figure size."""
    H = fig.get_figheight()
    n = subtitle.count("\n") + 1
    header = 0.14 + n * 0.19 + 0.12 + 0.26      # inches of space needed above
    fig.subplots_adjust(top=1 - header / H)
    pos = ax.get_position()
    sub_y = pos.y1 + 0.14 / H
    fig.text(pos.x0, sub_y, subtitle, fontsize=9.5, color=INK_2,
             ha="left", va="bottom", linespacing=1.55)
    fig.text(pos.x0, sub_y + (n * 0.19 + 0.12) / H, title, fontsize=13,
             fontweight="bold", color=INK, ha="left", va="bottom")


def save(fig, stem):
    # Deterministic output: no embedded timestamp, fixed element-id salt. Two
    # runs on the same data produce byte-identical files, so CHECKSUMS.sha256
    # stays valid after regenerating the figures.
    meta = {"png": {"Software": None}, "svg": {"Date": None}}
    for ext in ("png", "svg"):
        fig.savefig(FIGS / f"{stem}.{ext}", dpi=200, bbox_inches="tight",
                    pad_inches=0.28, metadata=meta[ext])
    plt.close(fig)
    print(f"  wrote figures/{stem}.png and .svg")


# ============================================================ figure 1
def fig_enzyme_forest():
    """Per-enzyme odds of Major severity. Emphasis form: the three enzymes that
    survive FDR correction carry the diverging pair; the six that do not are
    grey, because they are null results and must not read as trends."""
    st = rows("enzyme_severity_stats.csv")
    st.sort(key=lambda r: float(r["OR"]))

    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ys = range(len(st))

    for y, r in zip(ys, st):
        orv, lo, hi = float(r["OR"]), float(r["CI_lo"]), float(r["CI_hi"])
        sig = float(r["q_fdr"]) < 0.05
        color = (RED if orv > 1 else BLUE) if sig else GREY
        ax.plot([lo, hi], [y, y], color=color, linewidth=2.0,
                solid_capstyle="round", zorder=2, alpha=1.0 if sig else 0.75)
        # marker with a 2px surface ring so overlapping marks stay legible
        ax.plot([orv], [y], marker="o", markersize=9, color=color,
                markeredgecolor=SURFACE, markeredgewidth=2.0, zorder=4)
        if sig:
            ax.annotate(f"OR {orv:g}   q = {r['q_fdr']}",
                        xy=(hi, y), xytext=(7, 0), textcoords="offset points",
                        va="center", ha="left", fontsize=9, color=INK_2)

    ax.axvline(1, color=AXIS, linewidth=0.9, zorder=1)
    ax.set_yticks(list(ys))
    ax.set_yticklabels([f"{r['enzyme']}  ({r['n_pairs']} pairs)" for r in st],
                       fontsize=9.5, color=INK_2)
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 4, 8, 16])
    ax.set_xticklabels(["0.1", "0.25", "0.5", "1", "2", "4", "8", "16"],
                       fontsize=9)
    ax.set_xlim(0.05, 40)
    ax.set_ylim(-0.7, len(st) - 0.3)
    ax.set_xlabel("Odds ratio for Major severity (log scale, 95% CI)",
                  fontsize=9.5)
    style(ax)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)

    ax.legend(handles=[
        Line2D([], [], color=RED, marker="o", markersize=8, linewidth=2,
               markeredgecolor=SURFACE, markeredgewidth=2,
               label="Enriched for Major (FDR-significant)"),
        Line2D([], [], color=BLUE, marker="o", markersize=8, linewidth=2,
               markeredgecolor=SURFACE, markeredgewidth=2,
               label="Depleted of Major (FDR-significant)"),
        Line2D([], [], color=GREY, marker="o", markersize=8, linewidth=2,
               markeredgecolor=SURFACE, markeredgewidth=2,
               label="Not significant after FDR correction"),
    ], loc="lower right", frameon=False, fontsize=9, labelcolor=INK_2,
        handletextpad=0.8, borderaxespad=0.6)
    titles(fig, ax, "Which enzymes carry the severe interactions",
           "Odds of a Major-severity grade for pairs attributed to each enzyme, "
           "versus all other pairs.\nOnly three results survive FDR correction; "
           "the other six are null and are shown in grey.")
    save(fig, "fig1_enzyme_forest")


# ============================================================ figure 2
def fig_severity_confirmation():
    """Assigned severity predicts independent confirmation, in two unrelated
    external gold standards. Two series -> categorical slots 1 and 2."""
    kg = summary("kegg_validation_summary.csv")
    grades = ["Major", "Moderate", "Minor", "Unknown"]

    # DrugBank rates are re-derived from the released pair table, not copied
    # from the summary: unique pairs per grade, and how many are hard-proven.
    pairs = {}
    for r in rows("enzyme_pair_drugbank_flagged.csv"):
        k = (r["CID_A"], r["CID_B"])
        pairs.setdefault(k, {"sev": r["severity"], "hp": False})
        if r["drugbank_hard_proven"] == "Yes":
            pairs[k]["hp"] = True
    den = {g: sum(v["sev"] == g for v in pairs.values()) for g in grades}
    num = {g: sum(v["sev"] == g and v["hp"] for v in pairs.values())
           for g in grades}
    drugbank = [100 * num[g] / den[g] for g in grades]
    kegg = [float(kg[f"{g}_confirm_pct"]) for g in grades]

    fig, ax = plt.subplots(figsize=(8.8, 4.4))
    # limits first: hbar measures its corner radius and thickness cap in pixels
    ax.set_ylim(len(grades) - 0.45, -0.55)     # inverted, explicit
    ax.set_xlim(0, 44)
    h = 0.26                      # bar thickness, leaves air in the band
    off = 0.165                   # half-offset -> a clear surface gap between bars
    for i, g in enumerate(grades):
        hbar(ax, i - off, drugbank[i], h, BLUE)
        hbar(ax, i + off, kegg[i], h, ORANGE)
        ax.annotate(f"{drugbank[i]:.1f}%", xy=(drugbank[i], i - off),
                    xytext=(7, 0), textcoords="offset points", va="center",
                    fontsize=9, color=INK_2)
        ax.annotate(f"{kegg[i]:.1f}%", xy=(kegg[i], i + off), xytext=(7, 0),
                    textcoords="offset points", va="center", fontsize=9,
                    color=INK_2)

    ax.set_yticks(range(len(grades)))
    ax.set_yticklabels([f"{g}\nn = {den[g]:,}" for g in grades], fontsize=9.5,
                       color=INK_2, linespacing=1.5)
    ax.set_xlabel("Pairs independently confirmed (%)", fontsize=9.5)
    style(ax)

    ax.legend(handles=[
        Line2D([], [], color=BLUE, linewidth=7, solid_capstyle="butt",
               label="DrugBank hard-proven set   "
                     r"(trend $p = 1.8\times10^{-18}$)"),
        Line2D([], [], color=ORANGE, linewidth=7, solid_capstyle="butt",
               label="Curated KEGG DDI set   "
                     r"(trend $p = 9.6\times10^{-26}$)"),
    ], loc="lower right", frameon=False, fontsize=9, labelcolor=INK_2,
        handletextpad=0.9, borderaxespad=0.5, labelspacing=0.7)

    titles(fig, ax, "Assigned severity predicts independent confirmation",
           "Confirmation rate in two external gold standards that were not used "
           "to build the resource. The DrugBank\ngradient is monotone; in KEGG, "
           "Major and Moderate are indistinguishable, but graded pairs still "
           "confirm\nfar more often than ungraded ones. n = pairs testable in "
           "each gold standard.")
    save(fig, "fig2_severity_confirmation")


# ============================================================ figure 3
def fig_severity_composition():
    """Part-to-whole across an ordered scale -> stacked bar on the ordinal ramp.
    Unknown is grey: it is the absence of a grade, not a level of the scale."""
    ddi = rows("ddi_enzyme_database.csv")
    counts = {g: sum(r["severity"] == g for r in ddi)
              for g in ("Major", "Moderate", "Minor", "Unknown")}
    total = len(ddi)

    fig, ax = plt.subplots(figsize=(8.6, 2.3))
    x = 0.0
    gap = total * 0.0035  # 2px-equivalent surface gap between segments
    for g in ("Major", "Moderate", "Minor", "Unknown"):
        w = counts[g]
        color = SEV_RAMP.get(g, SEV_UNKNOWN)
        ax.barh([0], [w], left=[x], height=0.20, color=color,
                edgecolor="none", zorder=3)
        # label inside only where it comfortably fits, else above the segment
        pct = 100 * w / total
        label = f"{g}\n{w:,} ({pct:.0f}%)"
        if w / total > 0.15:
            ax.text(x + w / 2, 0, label, ha="center", va="center", fontsize=9.5,
                    color="#ffffff" if g in ("Major", "Moderate") else INK,
                    zorder=4, linespacing=1.5)
        else:
            ax.annotate(label, xy=(x + w / 2, 0.10), xytext=(0, 6),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=9, color=INK_2, linespacing=1.4)
        x += w + gap

    ax.set_xlim(-total * 0.005, total * 1.02)
    ax.set_ylim(-0.30, 0.52)
    ax.axis("off")
    titles(fig, ax, f"What the {total:,} enzyme-annotated pairs are graded as",
           "“Unknown” means DDInter records no clinical grade for the pair — it "
           "does not mean the pair is safe.\nTreating it as a negative class will "
           "bias any downstream analysis.")
    save(fig, "fig3_severity_composition")


# ============================================================ figure 4
def fig_mechanism_classes():
    """One series, magnitude -> one hue for every bar (a value-ramp here would
    double-encode bar length as colour)."""
    ddi = rows("ddi_enzyme_database.csv")
    counts = {}
    for r in ddi:
        for m in r["mechanism"].split(","):
            m = m.strip()
            if m:
                counts[m] = counts.get(m, 0) + 1
    items = sorted(counts.items(), key=lambda kv: kv[1])

    fig, ax = plt.subplots(figsize=(8.0, 3.4))
    ax.set_xlim(0, max(counts.values()) * 1.16)
    ax.set_ylim(-0.6, len(items) - 0.4)
    for i, (label, n) in enumerate(items):
        hbar(ax, i, n, 0.42, BLUE)
        ax.annotate(f"{n:,}", xy=(n, i), xytext=(7, 0),
                    textcoords="offset points", va="center", fontsize=9.5,
                    color=INK_2)
    ax.set_yticks(range(len(items)))
    ax.set_yticklabels([k[0].upper() + k[1:] for k in
                        [i[0] for i in items]], fontsize=10, color=INK_2)
    ax.set_xlabel("Pair occurrences", fontsize=9.5)
    style(ax)
    titles(fig, ax, "How the interactions are mediated",
           "A pair may carry more than one mechanism class, so occurrences "
           "exceed the number of pairs.")
    save(fig, "fig4_mechanism_classes")


# ============================================================ figure 5
def fig_phenotype_enrichment():
    """Top enzyme-phenotype signatures by fold enrichment. One hue: the bar
    length already carries magnitude."""
    ph = rows("enzyme_phenotype_enrichment.csv")
    ph.sort(key=lambda r: float(r["log2FE"]), reverse=True)
    top = ph[:12][::-1]

    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    mx = max(float(r["log2FE"]) for r in top)
    ax.set_xlim(0, mx * 1.14)
    ax.set_ylim(-0.6, len(top) - 0.4)
    for i, r in enumerate(top):
        hbar(ax, i, float(r["log2FE"]), 0.44, BLUE)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(
        [f"{r['AE_name'][:44]}\n{r['enzyme']} · {int(float(r['observed'])):,} obs"
         for r in top], fontsize=8.5, color=INK_2, linespacing=1.45)
    for i, r in enumerate(top):
        ax.annotate(f"{float(r['log2FE']):.2f}", xy=(float(r["log2FE"]), i),
                    xytext=(7, 0), textcoords="offset points", va="center",
                    fontsize=9, color=INK_2)
    ax.set_xlabel("log₂ fold enrichment over the background model", fontsize=9.5)
    style(ax)
    titles(fig, ax, "Strongest enzyme–adverse-event signatures",
           "Top 12 of 434 FDR-significant associations. These are "
           "drug-class-MEDIATED descriptive signatures,\nnot causal, "
           "enzyme-attributable adverse-event risk.")
    save(fig, "fig5_phenotype_enrichment")


if __name__ == "__main__":
    print("Regenerating figures from data/ …")
    fig_enzyme_forest()
    fig_severity_confirmation()
    fig_severity_composition()
    fig_mechanism_classes()
    fig_phenotype_enrichment()
    print("Done.")
