#!/usr/bin/env python3
"""Figure 1 of the refreshment-designs paper, from sims/sim1_results.csv (Simulation 1, true h = 0.2 cells).

For each selection strength beta and retention rate: the mean over 200 replications of the RELAXED plug-in set
(trimmed support, slack 0.10; translucent bar with end ticks), the population identified set of Theorem 3 for the
same data-generating process (dark segment; a point when beta = 0), and the mean naive quantile-alignment
estimate (dot).  Run from this folder: python3 make_fig1_funnel.py  (writes fig1_funnel.png and fig1_funnel.pdf).
"""
import csv, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
rows = list(csv.DictReader(open(HERE / "../../sims/sim1_results.csv")))
rows = [r for r in rows if abs(float(r["h_true"]) - 0.2) < 1e-9]
BLUE, AMBER, INK, MUT = "#2b6cb0", "#b45309", "#1a202c", "#5a6472"
COL = {"0.6": AMBER, "0.8": BLUE}
OFF = {"0.6": -0.05, "0.8": 0.05}

fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=200)
ax.axhline(0.2, ls="--", lw=1.1, color=MUT, zorder=1)
ax.annotate("true $h$ = 0.2", xy=(-0.27, 0.208), fontsize=8.6, color=MUT, va="bottom", ha="left")
for r in rows:
    b = float(r["beta"]); p = r["p_target"]
    x = b + OFF[p]; c = COL[p]
    hL, hU = float(r["mean_hL"]), float(r["mean_hU"])
    ax.plot([x, x], [hL, hU], color=c, lw=6, alpha=0.30, solid_capstyle="butt", zorder=2)
    ax.plot([x - 0.015, x + 0.015], [hL, hL], color=c, lw=1.8, zorder=3)
    ax.plot([x - 0.015, x + 0.015], [hU, hU], color=c, lw=1.8, zorder=3)
    pL, pU = float(r["pop_hL"]), float(r["pop_hU"])
    if pU - pL > 1e-9:
        ax.plot([x, x], [pL, pU], color=INK, lw=2.2, solid_capstyle="butt", zorder=4)
    else:
        ax.scatter([x], [pL], s=95, marker="s", facecolors="none", edgecolors=INK, linewidths=1.6, zorder=6)
    ax.scatter([x], [0.2 + float(r["naive_bias"])], s=30, color=c, zorder=5, edgecolors="white", linewidths=0.9)
ax.annotate("retention 60%", xy=(0 + OFF["0.6"] - 0.03, 0.47), fontsize=8.6, color=AMBER, ha="right", va="bottom")
ax.annotate("retention 80%", xy=(0 + OFF["0.8"] + 0.03, 0.06), fontsize=8.6, color=BLUE, ha="left", va="bottom")
ax.annotate("bar: relaxed plug-in set (trim at 10% of the mode, slack 0.10)\n"
            "dark segment (open square: a single point): population identified set (Theorem 3)\n"
            "dot: naive quantile-alignment estimate",
            xy=(0.02, 0.985), xycoords="axes fraction", fontsize=7.6, color=INK, va="top")
ax.set_xticks([0, 0.3, 0.6])
ax.set_xticklabels(["$\\beta = 0$\n(outcome-independent\nattrition)", "$\\beta = 0.3$", "$\\beta = 0.6$"], fontsize=8.6, color=INK)
ax.set_xlim(-0.3, 0.8)
ax.set_ylim(-0.16, 0.95)
ax.set_ylabel("conditioning shift $h$", fontsize=9.5, color=INK)
ax.spines[["top", "right"]].set_visible(False)
ax.spines[["left", "bottom"]].set_color("#cbd5e0")
ax.tick_params(colors=MUT, labelsize=8.6)
ax.grid(axis="y", color="#edf1f5", lw=0.8, zorder=0)
ax.set_title("Relaxed plug-in sets and population identified sets", fontsize=11, color=INK, loc="left", pad=12)
ax.annotate("Simulation 1, true $h$ = 0.2; means over 200 replications ($n$ = 4,000 + 2,000 refreshment)",
            xy=(0, 1.02), xycoords="axes fraction", fontsize=7.8, color=MUT)
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(HERE / f"fig1_funnel.{ext}", bbox_inches="tight", metadata={"CreationDate": None} if ext == "pdf" else None)
print("written")
