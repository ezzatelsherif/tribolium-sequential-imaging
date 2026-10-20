#!/usr/bin/env python3
"""Compare completed historical and confirmed-correction peak-order runs."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('historical', type=Path)
    parser.add_argument('corrected', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tables = {
        'temporal_peaks': 'work/figure7_ver05_analysis/gap_peaks/temporal_peak_references.csv',
        'pairrule_temporal_peaks': 'work/figure7_ver05_analysis/pairrule/temporal_peak_references.csv',
        'temporal_contrasts': 'work/figure7_ver05_analysis/gap_peaks/temporal_pairwise_contrasts.csv',
        'pairrule_temporal_contrasts': 'work/figure7_ver05_analysis/pairrule/all_temporal_pairwise_contrasts.csv',
        'comparison_inventory': 'S1_Data_Peak_Order_Comparisons.csv',
        'figure7_contrasts': 'Fig7_peak_time_differences_ver08.csv',
        'figure7_markers': 'Fig7_peak_markers_ver08.csv',
    }
    summary, changes = {}, []
    for name, relative in tables.items():
        a, b = [pd.read_csv(root / relative) for root in (args.historical, args.corrected)]
        if list(a.columns) != list(b.columns) or len(a) != len(b):
            raise ValueError(f'Incompatible table structure: {name}')
        columns = {}
        for column in a.columns:
            if column in {'source', 'source_file'}:
                continue
            x, y = a[column], b[column]
            if pd.api.types.is_numeric_dtype(x) and not pd.api.types.is_bool_dtype(x):
                same = np.isclose(x, y, rtol=0, atol=1e-12, equal_nan=True)
                finite = np.isfinite(x) & np.isfinite(y)
                delta = float(np.max(np.abs(y[finite] - x[finite]))) if finite.any() else 0.
            else:
                same = x.fillna('<missing>').eq(y.fillna('<missing>')).to_numpy()
                delta = None
            if not same.all():
                columns[column] = {'changed_rows': int((~same).sum()), 'max_absolute_difference': delta}
                for index in np.flatnonzero(~same):
                    key = {k: str(a.iloc[index][k]) for k in
                           ['pair_id', 'panel', 'cohort', 'event_id', 'stage', 'embryo', 'event_a', 'event_b']
                           if k in a.columns}
                    changes.append({'table': name, 'row_number': int(index)+1,
                                    'row_key': json.dumps(key), 'column': column,
                                    'historical': x.iloc[index], 'corrected': y.iloc[index]})
        summary[name] = {'rows': len(a), 'changed_columns': columns}
    profiles = 'work/figure7_ver05_analysis/gap_peaks/temporal_profiles.csv'
    a, b = [pd.read_csv(root / profiles) for root in (args.historical, args.corrected)]
    profile_rows = []
    for (cohort, gene), group in a.groupby(['cohort', 'gene'], sort=False):
        other = b.loc[group.index]
        for bound in [11., 15., 23.]:
            selected = group.substage_index.le(bound)
            profile_rows.append({'cohort': cohort, 'gene': gene, 'through_substage_index': bound,
                'max_smoothed_intensity_difference': float(np.max(np.abs(
                    group.loc[selected, 'smoothed_intensity'] - other.loc[selected, 'smoothed_intensity']))),
                'max_normalized_intensity_difference': float(np.max(np.abs(
                    group.loc[selected, 'normalized_intensity'] - other.loc[selected, 'normalized_intensity'])))})
    pd.DataFrame(profile_rows).to_csv(args.output_dir / 'profile_differences.csv', index=False)
    pd.DataFrame(changes).to_csv(args.output_dir / 'table_differences.csv', index=False)
    (args.output_dir / 'comparison_summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
