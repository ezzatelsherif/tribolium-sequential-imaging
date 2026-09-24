"""Recreate the paper's numerical profile panels and their underlying tables.

The spatial and temporal transformations follow Theophilus Frimpong's original
viewers. Figure layout, labels, and output handling are supplied here so that
the numerical panels can be generated without an interactive application.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from .data import (STAGES, GENE_FILES, COLORS, figure_config,
                   spatial_profiles, temporal_measurements,
                   summarize_temporal, comparison_temporal, output_directory)


def configured_spatial(config, entry, genes=None, blocks=False):
    """Read a selected embryo with the settings recorded for this run."""
    fraction = config["block_smoothing_fraction" if blocks else "spatial_smoothing_fraction"]
    return spatial_profiles(entry["path"], fraction, config["spatial_normalization"],
                            entry["genes"] if genes is None else genes,
                            config["spatial_smoothing_order"])


def save_figure(fig, directory, name):
    for suffix in ("png", "svg", "pdf"):
        fig.savefig(directory / f"{name}.{suffix}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_lines(ax, curves, genes, temporal=False):
    for gene in genes:
        x, y = curves[gene]
        ax.plot(x, y, color=COLORS[gene], lw=1.3,
                linestyle="--" if gene.endswith("int") else "-", label=gene)
    ax.set_ylim(-0.08, 1.08)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("Relative intensity")
    ax.set_xlabel("eve-en substage" if temporal else "AP position (0 anterior, 1 posterior)")


def block_intensity(x, y, gene):
    """Original viewer's contrast-enhanced block mapping on 800 positions.

    This mapping changes only the color display. It is not applied to the line
    plots or used to calculate measurements. The settings are fixed to those
    in the source viewer: alpha=0.3, contrast threshold=0.175, and gene-specific
    neighborhood half-widths of 75 (eve), 65 (runt/odd), or 50 positions.
    """
    dense_x = np.linspace(0, 1, 800)
    y = np.clip(np.interp(dense_x, x, np.clip(y, 0, 1)), 0, 1)
    window = {"eve": 75, "runt": 65, "odd": 65}.get(gene, 50)
    contrast = np.zeros_like(y)
    for index, value in enumerate(y):
        neighbors = np.r_[y[max(0, index-window):index],
                          y[index+1:min(len(y), index+window+1)]]
        minimum = neighbors.min()
        baseline = np.sort(neighbors)[:10].mean() if value > minimum + 0.0025 else minimum
        contrast[index] = abs(value - baseline)
    weight = np.where(y > 0.75, 1.0, np.minimum(contrast / 0.175, 1.0))
    alpha = np.where(y > 0.05, 0.9, 0.3) if gene == "cad" else 0.3
    intensity = np.where(y > 0.1, alpha*y + (1-alpha)*weight, y)
    return dense_x, np.clip(intensity, 0, 1)


def plot_blocks(ax, curves, genes, temporal=False, cutoff=24):
    rows = []
    for gene in genes:
        x, y = curves[gene]
        if temporal:
            x = (x - 1) / 23
        _, intensity = block_intensity(x, y, gene)
        color = np.asarray(to_rgb(COLORS[gene]))
        rows.append(1 - intensity[:, None] * (1 - color))
    ax.imshow(np.asarray(rows), aspect="auto", interpolation="nearest",
              extent=[1 if temporal else 0, 24 if temporal else 1, len(genes), 0])
    ax.set_yticks(np.arange(len(genes)) + 0.5, genes, fontsize=8)
    ax.set_xticks([])
    ax.tick_params(length=0)
    ax.set_xlim((1, cutoff) if temporal else (0, 1))


def comparison_figure(config, entries, temporal, output, figure_name, groups):
    """Create one four-stage comparison per gene group, including block displays."""
    dense_t = np.linspace(1, 24, 400)
    time_curves = {gene: (dense_t, PchipInterpolator(summary.stage_index,
                    summary.normalized_mean)(dense_t)) for gene, summary in temporal.items()}
    # Blocks in the original viewer use linear interpolation of stage means.
    time_blocks = {gene: (summary.stage_index.to_numpy(), summary.normalized_mean.to_numpy())
                   for gene, summary in temporal.items()}
    for panel, genes in groups.items():
        fig = plt.figure(figsize=(11, 12), layout="constrained")
        grid = fig.add_gridspec(8, 2, height_ratios=[3, 1] * 4)
        for row, stage in enumerate(config["comparison_stages"]):
            entry = next(item for item in entries if item["stage"] == stage)
            curves = configured_spatial(config, entry, genes)
            blocks = configured_spatial(config, entry, genes, blocks=True)
            cutoff = STAGES.index(stage) + 1
            spatial_ax = fig.add_subplot(grid[2*row, 0])
            temporal_ax = fig.add_subplot(grid[2*row, 1])
            plot_lines(spatial_ax, curves, genes)
            plot_lines(temporal_ax, time_curves, genes, temporal=True)
            spatial_ax.set_xlim(0, 1)
            temporal_ax.set_xlim(1, cutoff)
            temporal_ax.set_xticks(np.arange(1, cutoff+1), STAGES[:cutoff], rotation=60, fontsize=7)
            spatial_ax.set_title(f"Stage {stage}, {entry['embryo']}", fontsize=10)
            temporal_ax.set_title(f"Posterior history through {stage}", fontsize=10)
            plot_blocks(fig.add_subplot(grid[2*row+1, 0]), blocks, genes)
            plot_blocks(fig.add_subplot(grid[2*row+1, 1]), time_blocks, genes, True, cutoff)
        spatial_ax.legend(loc="lower left", fontsize=8, ncol=len(genes))
        fig.suptitle(f"{figure_name.replace('_', ' ')}{panel}: temporal-to-spatial comparison")
        save_figure(fig, output, f"{figure_name}{panel}")


def reproduce_profiles(output):
    output = output_directory(output) / "profiles"
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "svg.fonttype": "none", "pdf.fonttype": 42})
    config = figure_config()
    measurements = temporal_measurements()
    measurements.to_csv(output / "temporal_measurements.csv", index=False)
    summaries = []
    for (cohort, channel), _ in measurements.groupby(["cohort", "channel"]):
        for quantity in ("corrected", "workbook_display"):
            summaries.append(summarize_temporal(measurements, cohort, channel, quantity))
    pd.concat(summaries, ignore_index=True).to_csv(output / "temporal_stage_summary.csv", index=False)
    temporal = comparison_temporal(measurements, config["eve_reference"], config["temporal_plot_values"])
    pd.concat([table.assign(plotted_gene=gene) for gene, table in temporal.items()],
              ignore_index=True).to_csv(output / "comparison_temporal_profiles.csv", index=False)

    # Figure 5: each gene and its own co-measured eve reference. run/odd share eve.
    fig, axes = plt.subplots(7, 1, figsize=(11, 15), layout="constrained", sharex=True)
    dense = np.linspace(1, 24, 400)
    for ax, gene in zip(axes, GENE_FILES):
        for channel, plotted, color in (("mrna", gene, COLORS[gene]), ("eve", "eve", COLORS["eve"])):
            table = summarize_temporal(measurements, gene, channel, config["temporal_plot_values"])
            y = PchipInterpolator(table.stage_index, table.normalized_mean)(dense)
            sem = PchipInterpolator(table.stage_index, table.normalized_sem)(dense)
            ax.plot(dense, y, color=color, label=plotted)
            ax.fill_between(dense, y-sem, y+sem, color=color, alpha=0.2)
            ax.scatter(table.stage_index, table.normalized_mean, color=color, s=8)
        ax.set_ylabel("Relative intensity")
        ax.set_title(f"{gene}: posterior measurements", fontsize=10)
        ax.legend(loc="upper right", ncol=2, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xticks(range(1, 25), STAGES, rotation=60)
    axes[-1].set_xlabel("eve-en substage")
    save_figure(fig, output, "Figure_5_temporal_profiles")

    spatial_rows = []
    for entry in config["spatial_files"]:
        curves = configured_spatial(config, entry)
        for gene, (x, y) in curves.items():
            spatial_rows.append(pd.DataFrame({"source_file": entry["path"], "stage": entry["stage"],
                "embryo": entry["embryo"], "gene": gene, "ap_position": x, "normalized_smoothed": y}))
    pd.concat(spatial_rows, ignore_index=True).to_csv(output / "spatial_profiles.csv", index=False)

    for replica, name in ((1, "Figure_6"), (2, "Figure_S1"), (3, "Figure_S2")):
        entries = [item for item in config["spatial_files"]
                   if item["assay"] == "exonic" and item["embryo"] == f"EMB_{replica}"]
        groups = config["gene_groups"] if replica == 1 else config["supplement_gene_groups"]
        comparison_figure(config, entries, temporal, output, name, groups)

    intronic = [item for item in config["spatial_files"] if item["assay"] == "intron_exon"]
    genes = ["kr", "krint", "mlpt", "mlptint", "cad"]
    for name, stages, replicates in (("Figure_4B", ["3.3", "4.2", "4.3"], [1]),
                                      ("Figure_S5", ["4.2", "4.3"], [1, 3])):
        fig, axes = plt.subplots(len(stages), len(replicates), figsize=(6*len(replicates), 3*len(stages)),
                                 squeeze=False, layout="constrained")
        for i, stage in enumerate(stages):
            for j, embryo in enumerate(replicates):
                entry = next(item for item in intronic if item["stage"] == stage
                             and item["embryo"] == f"EMB_{embryo}")
                plot_lines(axes[i, j], configured_spatial(config, entry, genes), genes)
                axes[i, j].set_title(f"Stage {stage}, replicate {j+1} ({entry['embryo']})")
                axes[i, j].set_xlim(0, 1)
        axes[0, 0].legend(ncol=3, fontsize=8)
        save_figure(fig, output, name)

    example = next(item for item in config["spatial_files"]
                   if "Figure 2C" in item["figures"])
    curves = configured_spatial(config, example)
    fig, ax = plt.subplots(figsize=(11, 4), layout="constrained")
    plot_lines(ax, curves, list(COLORS.keys())[:10])
    ax.legend(ncol=5, fontsize=8)
    ax.set_title("Figure 2C: spatial profiles from stage 4.1, EMB_1")
    save_figure(fig, output, "Figure_2C_profiles")

    example = next(item for item in config["spatial_files"]
                   if "Figure 3F" in item["figures"])
    groups = [["kr", "mlpt", "gt", "svb", "cad"], ["kr", "krint"],
              ["mlpt", "mlptint"], ["svb", "svbint"]]
    curves = configured_spatial(config, example)
    fig, axes = plt.subplots(4, 1, figsize=(10, 11), layout="constrained")
    for ax, genes in zip(axes, groups):
        plot_lines(ax, curves, genes)
        ax.set_xlim(0, 1)
        ax.legend(ncol=len(genes), fontsize=8)
    fig.suptitle("Figure 3F: spatial overlays from stage 4.2, EMB_1")
    save_figure(fig, output, "Figure_3F_spatial_overlays")
    (output / "settings.json").write_text(json.dumps(config, indent=2) + "\n")
    print(f"Wrote profiles, tables, and plotting settings to {output}")
