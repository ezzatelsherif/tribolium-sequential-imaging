"""Read the paper's extracted fluorescence measurements without modifying them.

Spatial files contain AP bins from one embryo. Temporal workbook rows contain
independent embryos at one assigned eve-en substage. Numbered specimen slots
are paired only within a workbook row, never across developmental stages.
"""
from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

ROOT = Path(__file__).resolve().parents[1]
STAGES = tuple(f"{cycle}.{phase}" for cycle in range(1, 9) for phase in (1, 2, 3))
GENE_FILES = {"hb": "HB.xlsx", "kr": "KR.xlsx", "mlpt": "MLPT.xlsx",
              "gt": "GT.xlsx", "svb": "SVB.xlsx", "runt": "run.xlsx",
              "odd": "odd.xlsx"}
COLORS = {"hb": "blue", "kr": "red", "mlpt": "green", "gt": "gold",
          "svb": "brown", "eve": "purple", "runt": "hotpink", "odd": "cyan",
          "en": "olive", "cad": "violet", "krint": "#ff8080",
          "mlptint": "#80c080", "svbint": "#c08080", "eve-en": "purple"}


def figure_config():
    """Return the recorded source-file assignments and plotting settings."""
    return json.loads((ROOT / "config/figures.json").read_text())


def normalize(values, method="minmax"):
    """Apply the recorded display normalization to a finite one-dimensional array."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all() or not len(values):
        raise ValueError("Expected a nonempty, finite, one-dimensional profile")
    baseline = values.min() if method == "minmax" else 0.0
    if method not in ("minmax", "maximum"):
        raise ValueError(f"Unknown normalization: {method}")
    scale = values.max() - baseline
    if scale <= 0:
        raise ValueError("The profile has no positive normalization range")
    return (values - baseline) / scale


def smooth(values, fraction=0.1, order=3):
    """Savitzky-Golay smoothing with Theo's original odd-window rule.

    The default window spans about 10% of the measured AP bins. Smoothing is
    performed after normalization and can produce small excursions beyond 0-1.
    Zero disables smoothing. Values are not clipped in line plots.
    """
    values = np.asarray(values, dtype=float)
    if not 0 <= fraction <= 1:
        raise ValueError("Smoothing fraction must be between 0 and 1")
    if fraction == 0 or len(values) < order + 2:
        return values.copy()
    n = len(values)
    window = int(max(order + 2, (n * fraction) // 2 * 2 + 1))
    window = min(window, n if n % 2 else n - 1)
    if window % 2 == 0:
        window += 1
    return savgol_filter(values, window_length=window, polyorder=order)


def spatial_profiles(relative_path, fraction=0.1, normalization="minmax", genes=None, order=3):
    """Read one embryo, retaining the correspondence between position and signal.

    Entirely empty trailing CSV rows are ignored. Missing individual values in
    requested channels are rejected. Other channels do not affect those curves.
    The original distance column is preserved in the source file. Its range is
    mapped to 0-1 for the spatial comparisons.
    """
    table = pd.read_csv(ROOT / relative_path).dropna(how="all")
    coordinate = "Distance_(microns)"
    if coordinate not in table:
        raise ValueError(f"Missing {coordinate}: {relative_path}")
    if genes is None:
        genes = [gene for gene in table.columns if gene != coordinate]
    values = table[[coordinate, *genes]].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"Incomplete spatial rows in {relative_path}")
    distance = table[coordinate].to_numpy(dtype=float)
    if not np.all(np.diff(distance) > 0):
        raise ValueError(f"AP coordinates must be strictly increasing: {relative_path}")
    x = normalize(distance)
    return {gene: (x, smooth(normalize(table[gene], normalization), fraction, order))
            for gene in genes}


def output_directory(path):
    """Resolve a run location while preventing writes inside source directories."""
    path = Path(path).resolve()
    reserved = [ROOT / name for name in ("analysis", "data", "seqimaging", "config", "validation")]
    if path == ROOT or any(path == source or source in path.parents for source in reserved):
        raise ValueError("Choose an output directory outside the repository's source directories")
    return path


def temporal_measurements():
    """Return all source slots, including incomplete measurements, in long form.

    A is the posterior ROI and B the anterior background ROI. Negative A-B
    values are retained. Completeness is evaluated separately for each channel;
    the paired intron/exon analysis applies its additional paired-slot filter.
    """
    rows = []
    for gene, filename in GENE_FILES.items():
        path = ROOT / "data/temporal" / filename
        table = pd.read_excel(path, sheet_name="Sheet1")
        displayed = pd.read_excel(path, sheet_name="BASE")
        displayed.columns = displayed.columns.str.strip()
        labels = [f"{value:.1f}" for value in table["Stage"]]
        if labels != list(STAGES):
            raise ValueError(f"Unexpected stage order in {filename}")
        for index, row in table.iterrows():
            for channel in ("eve", "mrna", "intr"):
                if f"{channel}1A" not in table:
                    continue
                for specimen in range(1, 6):
                    posterior = row[f"{channel}{specimen}A"]
                    background = row[f"{channel}{specimen}B"]
                    complete = bool(np.isfinite(posterior) and np.isfinite(background))
                    # Retain the cached values that the original plotting code read.
                    # A separate raw A-B column makes any source discrepancy visible.
                    plot_value = displayed.loc[index, f"{channel}{specimen}"]
                    if plot_value == 0:
                        plot_value = np.nan
                    rows.append({"cohort": gene, "stage": labels[index],
                                 "stage_index": index + 1, "specimen_slot": specimen,
                                 "channel": channel, "posterior": posterior,
                                 "anterior_background": background,
                                 "corrected": posterior - background,
                                 "workbook_display": plot_value,
                                 "complete": complete})
    return pd.DataFrame(rows)


def summarize_temporal(measurements, cohort, channel, quantity="corrected"):
    """Stage means and SEM, scaled by the full-series range of stage means."""
    selected = measurements.loc[(measurements.cohort == cohort)
                                & (measurements.channel == channel)
                                & np.isfinite(measurements[quantity])]
    grouped = selected.groupby("stage_index")[quantity]
    summary = grouped.agg(["count", "mean", "std"]).reindex(range(1, 25))
    if summary["count"].isna().any() or (summary["count"] < 2).any():
        raise ValueError(f"Insufficient measurements for {cohort}/{channel}")
    summary["sem"] = summary["std"] / np.sqrt(summary["count"])
    minimum, maximum = summary["mean"].min(), summary["mean"].max()
    if maximum <= minimum:
        raise ValueError(f"No temporal range for {cohort}/{channel}")
    summary["normalized_mean"] = (summary["mean"] - minimum) / (maximum - minimum)
    summary["normalized_sem"] = summary["sem"] / (maximum - minimum)
    summary["normalization_minimum"] = minimum
    summary["normalization_range"] = maximum - minimum
    summary["stage"] = STAGES
    summary["cohort"], summary["channel"] = cohort, channel
    summary["quantity"] = quantity
    return summary.reset_index()


def comparison_temporal(measurements, eve_reference, quantity="workbook_display"):
    """Prepare the temporal curves with an explicit eve-reference choice.

    'runt' selects the eve observations from the shared pair-rule cohort.
    'pooled' reproduces the original viewer's concatenation of the seven
    workbooks, including the identical eve entries in run.xlsx and odd.xlsx.
    Individual cohort summaries remain available separately in the output.
    """
    summaries = {gene: summarize_temporal(measurements, gene, "mrna", quantity)
                 for gene in GENE_FILES}
    if eve_reference == "pooled":
        pooled = measurements.loc[measurements.channel == "eve"].copy()
        pooled["cohort"] = "pooled"
        summaries["eve"] = summarize_temporal(pooled, "pooled", "eve", quantity)
    elif eve_reference in GENE_FILES:
        summaries["eve"] = summarize_temporal(measurements, eve_reference, "eve", quantity)
    else:
        raise ValueError(f"Unknown eve reference: {eve_reference}")
    return summaries
