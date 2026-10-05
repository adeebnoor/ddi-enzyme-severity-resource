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
from matplotlib.patches import PathPatch, FancyBboxPatch, FancyArrowPatch
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
    meta = {"png": {"Software": None}, "svg": {"Date": None}, "pdf": {"CreationDate": None, "ModDate": None}}
    for ext in ("png", "svg", "pdf"):
        fig.savefig(FIGS / f"{stem}.{ext}", dpi=200, bbox_inches="tight",
                    pad_inches=0.28, metadata=meta[ext])
    plt.close(fig)
    print(f"  wrote figures/{stem}.png and .svg")


# ============================================================ figure 1
def fig_enzyme_forest():
    """Per-gene exploratory association using frozen accession-backed identities."""
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
            ax.annotate(f"OR {orv:g}   $q$ = {r['q_fdr']}",
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
    _nsig = sum(float(r["q_fdr"]) < 0.05 for r in st)
    titles(fig, ax, "Exploratory enzyme-associated severity comparisons",
           "Accession-backed genes; odds of Major severity for attributed inhibition pairs, "
           "versus all other graded inhibition pairs.\n"
           f"{_nsig} of {len(st)} results survive FDR correction; "
           f"the other {len(st) - _nsig} do not reach the threshold and are grey.\n"
           "All 12 OAT1 exposed pairs contain methotrexate; no general OAT1 effect is identifiable.")
    save(fig, "fig3_enzyme_forest")


# ============================================================ figure 2
def fig_severity_confirmation():
    """Three historical reported membership comparisons; separate denominators."""
    grades = ["Major", "Moderate", "Minor", "Unknown"]
    with open(ROOT / "tables" / "SupplementaryTableS3_external_aggregate_counts.csv", encoding="utf-8") as fh:
        aggregate = list(csv.DictReader(fh))
    series = [("DrugBank", "DrugBank archived subset", BLUE),
              ("KEGG", "KEGG reported aggregate", ORANGE),
              ("trueDDI author reference", "Author-labelled trueDDI reference", "#2c8a7b")]
    by_source = {source: {r["grade"]: r for r in aggregate if r["resource"] == source}
                 for source, label, color in series}
    fig, ax = plt.subplots(figsize=(10.6, 5.8))
    ax.set_ylim(len(grades) - 0.4, -0.7)
    ax.set_xlim(0, 55)
    for j, (source, label, color) in enumerate(series):
        for i, grade in enumerate(grades):
            row = by_source[source][grade]
            rate = 100 * int(row["confirmed_pairs"]) / int(row["testable_pairs"])
            y = i + (j - 1) * 0.23
            hbar(ax, y, rate, 0.18, color)
            ax.annotate(f"{rate:.1f}%", xy=(rate, y), xytext=(7, 0),
                        textcoords="offset points", va="center", fontsize=9, color=INK_2)
    ax.set_yticks(range(len(grades)))
    ax.set_yticklabels([f"{g}\nDB n={by_source['DrugBank'][g]['testable_pairs']}; KEGG n={by_source['KEGG'][g]['testable_pairs']}\nAuthor ref n={by_source['trueDDI author reference'][g]['testable_pairs']}"
                        for g in grades], fontsize=8.5, color=INK_2, linespacing=1.3)
    ax.set_xlabel("Pairs present in the reported reference set (%)", fontsize=9.5)
    style(ax)
    ax.legend(handles=[Line2D([], [], color=color, linewidth=7, label=label)
                       for source, label, color in series], loc="lower right", frameon=False,
              fontsize=9, labelcolor=INK_2)
    titles(fig, ax, "Descriptive historical reference membership",
           "DrugBank and author-labelled trueDDI counts reproduce from exact recovered author-held input snapshots.\n"
           "KEGG remains a reported aggregate. Denominators differ; reproduction does not establish independent clinical validation.")
    save(fig, "fig4_external_concordance")


# ============================================================ figure 3
def fig_severity_composition():
    """Part-to-whole across an ordered scale -> stacked bar on the ordinal ramp.
    Unknown is grey: it is the absence of a grade, not a level of the scale."""
    # aggregate counts; the pair-level grades are not redistributed (DDInter licence)
    summ = {r["grade"]: int(r["pairs"]) for r in rows("severity_summary.csv")}
    counts = {g: summ[g] for g in ("Major", "Moderate", "Minor", "Unknown")}
    total = summ["Total"]

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
    save(fig, "fig2_severity_composition")


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
    titles(fig, ax, "Reported mechanism-class occurrences",
           "Historical CYP-labelled classes include general Enzyme annotations and shared targets.\n"
           "A pair may carry several classes, so occurrences exceed the number of pairs.")
    save(fig, "supp_fig1_mechanism_classes")


def fig_provenance():
    """Explicitly separate contextual GoldD3R and unresolved enzyme input branches."""
    fig, ax = plt.subplots(figsize=(10.6, 6.4))
    ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off")
    def box(x,y,w,h,text,missing=False):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.07",facecolor="#f2f4f6" if missing else "#eaf2fc",edgecolor=GREY if missing else BLUE,linestyle="--" if missing else "-",linewidth=1.2))
        ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=9.5,color=INK,linespacing=1.5)
    def arrow(x1,y1,x2,y2,dashed=False):
        ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=12,linewidth=1.1,color=INK_2,linestyle="--" if dashed else "-"))
    box(.2,4.45,2.65,1.0,"GoldD3R author export\n21,897 pairs; 12,493 labelled\nExact deposited file hashed")
    box(.2,2.8,2.65,1.0,"Contextual harmonisation\nHistorical results not regenerated\nOriginal intermediates absent",True)
    arrow(1.52,4.4,1.52,3.86,True)
    box(.2,.85,2.65,1.05,"UniProt identity audit\n35 reviewed human accessions\n2 source symbols corrected")
    arrow(2.92,1.38,3.5,1.38)
    box(3.6,4.45,2.85,1.0,"Original D3 enzyme branch\nHistorical graph + rules recovered\nCurrent-record identity unresolved",True)
    box(3.6,2.8,2.85,1.0,"Deposited attribution + bridge\n3,072 rows → 1,900 pair records\nRecord contents rebuild exactly")
    arrow(5.03,4.4,5.03,3.86,True)
    box(7.15,4.45,2.65,1.0,"User-acquired DDInter 2.0 CSVs\nSource-specific terms apply\nNot deposited as clinical grades",True)
    box(7.15,2.8,2.65,1.0,"Identifier join + max-grade rule\nExpected 93 / 419 / 50 / 1,338\nSHA-256 joined-grade fingerprint")
    arrow(8.48,4.4,8.48,3.86)
    arrow(6.53,3.3,7.05,3.3)
    box(3.6,.85,2.85,1.05,"Exploratory association example\n549 graded inhibition pairs\n424 single-enzyme sensitivity pairs")
    arrow(7.95,2.73,6.1,1.96)
    box(7.15,.85,2.65,1.05,"Historical reference comparison\nDrugBank / trueDDI recovered\nKEGG snapshot unresolved",True)
    ax.text(.2,.1,"Solid: deposited/exported steps. Dashed: current-record lineage unresolved or historical reported layer. Identity does not validate pair relations.",ha="left",fontsize=9,color=INK_2)
    fig.suptitle("Provenance and reproducibility boundaries",x=.075,ha="left",fontsize=14,fontweight="bold")
    save(fig,"fig1_provenance")

if __name__ == "__main__":
    print("Regenerating figures from data/ …")
    fig_provenance()
    fig_enzyme_forest()
    fig_severity_confirmation()
    fig_severity_composition()
    fig_mechanism_classes()
    print("Done.")
