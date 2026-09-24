"""Check data provenance and compare regenerated curves with manuscript traces."""
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from .data import (ROOT, figure_config, spatial_profiles, temporal_measurements,
                   comparison_temporal, normalize)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify():
    manifest = pd.read_csv(ROOT / "data/manifest.csv")
    require(not manifest.path.duplicated().any(), "Duplicate manifest paths")
    found = {str(p.relative_to(ROOT)) for suffix in ("*.csv", "*.xlsx")
             for p in (ROOT / "data").rglob(suffix)} - {"data/manifest.csv"}
    require(found == set(manifest.path), "Unlisted or missing source data")
    for row in manifest.itertuples():
        contents = (ROOT / row.path).read_bytes()
        require(len(contents) == row.bytes, f"File size changed: {row.path}")
        require(hashlib.sha256(contents).hexdigest() == row.sha256,
                f"Source checksum changed: {row.path}")
    config = figure_config()
    assignments = config["spatial_files"]
    assigned = [entry["path"] for entry in assignments]
    spatial_files = {p for p in found if p.startswith("data/spatial/")}
    require(len(assigned) == len(set(assigned)) == 17, "Expected 17 spatial source files")
    require(set(assigned) == spatial_files, "Spatial files and figure assignments differ")
    require(all(entry["figures"] for entry in assignments), "Unassigned source data")
    for name in ("KR.xlsx", "MLPT.xlsx", "SVB.xlsx"):
        original = ROOT / "analysis/intron_exon/empirical/temporal" / name
        require(original.read_bytes() == (ROOT / "data/temporal" / name).read_bytes(),
                f"The model and profile source workbooks differ: {name}")
    print(f"Verified {len(manifest)} source checksums and all 17 spatial figure assignments.")

    measurements = temporal_measurements()
    complete = measurements.loc[measurements.complete]
    require(np.allclose(complete.corrected, complete.workbook_display, atol=1e-9),
            "Cached workbook values differ from available posterior-minus-background measurements")
    exception = measurements.loc[~measurements.complete & measurements.workbook_display.notna()]
    require(len(exception) == 1, "The recorded missing-background exception changed")
    row = exception.iloc[0]
    require((row.cohort, row.stage, row.channel, row.specimen_slot, row.posterior,
             row.workbook_display) == ("kr", "8.1", "eve", 5, 1380, 1380),
            "Unexpected incomplete measurement used by the original viewer")
    require(np.isnan(row.anterior_background), "Expected a missing anterior background")
    for gene, count in (("kr", 116), ("mlpt", 118), ("svb", 116)):
        pairs = measurements.loc[(measurements.cohort == gene)
                                  & measurements.channel.isin(["mrna", "intr"])].pivot(
            index=["stage", "specimen_slot"], columns="channel", values="corrected")
        require(len(pairs.dropna()) == count, f"Unexpected paired specimen count for {gene}")
    eve = measurements.loc[measurements.channel == "eve"]
    run = eve.loc[eve.cohort == "runt", ["posterior", "anterior_background"]].to_numpy()
    odd = eve.loc[eve.cohort == "odd", ["posterior", "anterior_background"]].to_numpy()
    require(np.array_equal(run, odd, equal_nan=True), "The shared run/odd eve measurements differ")
    print("Verified raw corrections and paired specimen counts (Kr 116, mlpt 118, svb 116).")
    print("Recorded source exception: KR.xlsx, stage 8.1, eve slot 5 has no anterior background;")
    print("  the original display uses cached 1380; raw A-B remains missing. See docs/data_notes.md.")

    reference = json.loads((ROOT / "validation/figure_traces.json").read_text())
    x = np.linspace(0, 1, reference["points"])
    temporal = comparison_temporal(measurements, config["eve_reference"], config["temporal_plot_values"])
    spatial, errors = {}, []
    for trace in reference["curves"]:
        if trace["kind"] == "spatial":
            source = trace["source"]
            if source not in spatial:
                genes = next(entry["genes"] for entry in assignments if entry["path"] == source)
                spatial[source] = spatial_profiles(source, config["spatial_smoothing_fraction"],
                                                    config["spatial_normalization"], genes,
                                                    config["spatial_smoothing_order"])
            xx, yy = spatial[source][trace["gene"]]
        else:
            series = temporal[trace["gene"]]
            xx = np.linspace(0, 1, 400)
            yy = PchipInterpolator((series.stage_index-1)/23, series.normalized_mean)(xx)
        predicted = np.interp(x, xx, normalize(yy))
        error = float(np.sqrt(np.mean((predicted - trace["y"])**2)))
        require(error < reference["rmse_tolerance"],
                f"Figure curve changed: {trace['figure']}, {trace['gene']}, RMSE={error:.6f}")
        errors.append(error)
    require(set(spatial) == spatial_files, "Not every spatial embryo has a manuscript trace")
    print(f"Matched {len(errors)} reference curve shapes; maximum RMSE {max(errors):.6f}.")
