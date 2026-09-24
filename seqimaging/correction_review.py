"""Compare the historical eve reference with the confirmed data corrections.

Run from the repository root: python -m seqimaging.correction_review
Outputs are review tables and plots, separate from the manuscript panels.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .data import ROOT, STAGES, temporal_measurements, comparison_temporal


def review():
    source = temporal_measurements(apply_corrections=False)
    corrected = temporal_measurements()
    scenarios = {
        "historical": comparison_temporal(source, "pooled", "workbook_display")["eve"],
        "background_corrected": comparison_temporal(corrected, "pooled")["eve"],
        "pooled_unique": comparison_temporal(corrected, "pooled_unique")["eve"],
        "run_cohort": comparison_temporal(corrected, "runt")["eve"],
    }
    output = ROOT / "outputs/correction_review"
    output.mkdir(parents=True, exist_ok=True)
    pd.concat([table.assign(scenario=name) for name, table in scenarios.items()]).to_csv(
        output / "eve_reference_comparison.csv", index=False)
    baseline = scenarios["historical"].normalized_mean.to_numpy()
    metrics = {}
    for name, table in scenarios.items():
        y = table.normalized_mean.to_numpy()
        difference = y - baseline
        peaks = [STAGES[i] for i in range(1, 23) if y[i] > y[i-1] and y[i] > y[i+1]]
        metrics[name] = {
            "maximum_absolute_change": float(abs(difference).max()),
            "stage_of_maximum_change": STAGES[abs(difference).argmax()],
            "root_mean_square_change": float(np.sqrt(np.mean(difference**2))),
            "peak_substages": peaks,
            "minimum_stage_count": int(table["count"].min()),
            "maximum_stage_count": int(table["count"].max()),
        }
    (output / "comparison_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    plt.rcParams.update({"font.size": 10, "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.5), sharex=True, layout="constrained")
    labels = {"historical": "Historical seven-workbook pool",
              "pooled_unique": "Corrected pool, shared cohort counted once",
              "run_cohort": "Run cohort only"}
    colors = {"historical": "#555555", "pooled_unique": "#7b3294", "run_cohort": "#008837"}
    for name, label in labels.items():
        y = scenarios[name].normalized_mean.to_numpy()
        axes[0].plot(range(1, 25), y, marker="o", markersize=3,
                     color=colors[name], label=label, linestyle="--" if name == "historical" else "-")
        if name != "historical":
            axes[1].plot(range(1, 25), y-baseline, marker="o", markersize=3,
                         color=colors[name], label=label)
    axes[0].set_ylabel("Normalized eve stage mean")
    axes[0].legend(fontsize=9, loc="lower left", bbox_to_anchor=(0, 1.02), frameon=False)
    axes[1].axhline(0, color="#777777", linewidth=0.8)
    axes[1].set_ylabel("Change from historical curve")
    axes[1].set_xticks(range(1, 25), STAGES, rotation=60)
    axes[1].set_xlabel("eve-en substage")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(output / f"eve_reference_comparison.{suffix}", dpi=180)
    plt.close(fig)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    review()
