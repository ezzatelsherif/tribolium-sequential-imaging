"""Reconstruct the exact observed and jointly resampled stage means.

Uses the original load_sources and boot_cohorts functions unchanged; avoids
producing superseded raw half-height analyses and their large caches.
"""
import argparse
import json
import numpy as np
import pandas as pd
from temporal_direct import ROOT, load_sources, boot_cohorts


def main(count, apply_corrections=False):
    cohorts, audit = load_sources()
    corrections = []
    if apply_corrections:
        table = pd.read_csv(ROOT / 'data/temporal_corrections.csv', dtype={'stage': str})
        labels = [f'{cycle}.{stage}' for cycle in range(1, 9) for stage in range(1, 4)]
        for row in table.itertuples():
            cohort = row.cohort.upper() if row.cohort not in ('run', 'odd', 'pairrule') else row.cohort
            channel = row.channel
            stage, slot = labels.index(row.stage), int(row.specimen_slot)-1
            workbook = pd.read_excel(ROOT / f'data/Temporal/{cohort}.xlsx', sheet_name='Sheet1')
            prefix = channel if channel == 'eve' else 'mrna'
            source = workbook.iloc[stage]
            if source[f'{prefix}{slot+1}A'] != row.posterior or np.isfinite(source[f'{prefix}{slot+1}B']):
                raise ValueError(f'Correction source does not match: {row.correction_id}')
            if np.isfinite(cohorts[cohort][channel][stage, slot]):
                raise ValueError(f'Correction target was not missing: {row.correction_id}')
            value = float(row.posterior-row.anterior_background)
            cohorts[cohort][channel][stage, slot] = value
            corrections.append({'correction_id': row.correction_id, 'cohort': cohort,
                                'channel': channel, 'stage': row.stage,
                                'specimen_slot': slot+1, 'corrected_value': value})
    curves = boot_cohorts(cohorts, count)
    paired = {}
    for cohort, channels in cohorts.items():
        complete = np.stack([np.isfinite(a) for a in channels.values()]).all(0)
        paired[cohort] = {gene: np.where(complete, a, np.nan)
                          for gene, a in channels.items()}
    paired_curves = boot_cohorts(paired, count, '_paired')
    np.savez_compressed(ROOT / 'temporal_bootstrap_curves.npz',
        **{f'{cohort}_{gene}': values for cohort, channels in curves.items()
           for gene, values in channels.items()},
        **{f'paired_{cohort}_{gene}': values for cohort, channels in paired_curves.items()
           for gene, values in channels.items()})
    (ROOT / 'bootstrap_audit.json').write_text(json.dumps({
        'bootstrap_draws': count, 'seed': 20261019,
        'draw_zero': 'observed stage means',
        'independent_unit': 'embryo within assigned developmental stage',
        'paired_channels_resampled_together': True,
        'confirmed_corrections_applied': corrections,
        'source_audit': audit}, indent=2))
    print(f'Created observed means and {count} bootstrap draws from raw workbooks.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bootstraps', type=int, default=20000)
    parser.add_argument('--apply-confirmed-corrections', action='store_true')
    args = parser.parse_args()
    main(args.bootstraps, args.apply_confirmed_corrections)
