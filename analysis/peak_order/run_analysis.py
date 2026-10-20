#!/usr/bin/env python3
"""Reproduce the peak-order analysis from raw measurements in an isolated folder."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

PACKAGE = Path(__file__).resolve().parent


def validate_tables(work: Path) -> list[dict]:
    import pandas as pd
    results = []
    for expected in sorted((PACKAGE / 'reference').rglob('*.csv')):
        relative = expected.relative_to(PACKAGE / 'reference')
        actual = work / relative
        left, right = pd.read_csv(expected), pd.read_csv(actual)
        # Figure export provenance paths vary with the user's output directory.
        ignored = [c for c in left if c in {'source', 'source_file'}]
        pd.testing.assert_frame_equal(left.drop(columns=ignored), right.drop(columns=ignored),
                                      check_exact=False, rtol=1e-12, atol=1e-12)
        results.append({'table': str(relative), 'rows': len(left),
                        'columns_compared': len(left.columns)-len(ignored),
                        'ignored_path_columns': ignored, 'passed': True})
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=PACKAGE / 'results')
    parser.add_argument('--data-root', type=Path, default=PACKAGE.parents[1] / 'data',
                        help='Repository data directory (temporal/ and spatial/exonic/).')
    parser.add_argument('--apply-confirmed-corrections', action='store_true',
                        help='Sensitivity run: apply temporal_corrections.csv; does not replace the reviewed historical run.')
    parser.add_argument('--bootstraps', type=int, default=20000,
                        help='Use 20000 to reproduce the reported intervals; lower counts are smoke tests.')
    parser.add_argument('--skip-figures', action='store_true', help='Reproduce tables without rendering figures.')
    parser.add_argument('--no-reference-check', action='store_true',
                        help='Skip comparison with reviewed 20000-draw source tables.')
    args = parser.parse_args()
    if args.bootstraps < 2:
        parser.error('--bootstraps must be at least 2')
    data = args.data_root.expanduser().resolve()
    if not (data / 'temporal').is_dir() or not (data / 'spatial/exonic').is_dir():
        parser.error(f'Missing temporal/ or spatial/exonic/ under {data}')
    output = args.output_dir.expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        parser.error(f'Output directory is not empty: {output}. Choose a new directory.')
    output.mkdir(parents=True, exist_ok=True)
    work = output / 'work'
    shutil.copytree(PACKAGE / 'engine', work)
    destination = work / 'all_gene_direct_analysis/data'
    shutil.copytree(data / 'temporal', destination / 'Temporal')
    shutil.copytree(data / 'spatial/exonic', destination, dirs_exist_ok=True)
    if args.apply_confirmed_corrections:
        shutil.copyfile(data / 'temporal_corrections.csv', destination / 'temporal_corrections.csv')
    shutil.copyfile(PACKAGE / 'inputs/reconstructed_points.json',
                    work / 'all_gene_direct_analysis/reconstructed_points.json')
    shutil.copytree(PACKAGE / 'fonts', work / 'fonts')
    (work / 'figure7_ver04_analysis/kr_eve').mkdir(exist_ok=True)
    shutil.copyfile(work / 'figure7_ver08/plot_core.py', work / 'figure7_ver06/plot_core.py')
    logdir = output / 'logs'; logdir.mkdir()
    py = sys.executable
    steps = [
        ('01_bootstrap', 'all_gene_direct_analysis', [py, 'bootstrap_only.py', '--bootstraps', str(args.bootstraps)] +
         (['--apply-confirmed-corrections'] if args.apply_confirmed_corrections else [])),
        ('02_spatial_raw', 'all_gene_direct_analysis', [py, 'spatial_direct.py']),
        ('03_pairrule_smoothing', 'figure7_ver04_analysis/pairrule', [py, '-c',
         'import extract_pairrule as p; p.spatial_analysis(p.temporal_analysis()); p.annotate_neighbor_tails()']),
        ('04_gt_mlpt_smoothing', 'figure7_ver04_analysis/gt_mlpt', [py, 'reanalyze.py']),
        ('05_svb_kr_smoothing', 'figure7_ver04_analysis/svb_eve', [py, 'extract_svb_eve.py']),
        ('06_svb_sensitivity', 'figure7_ver04_analysis/svb_eve', [py, 'sensitivity_and_checks.py']),
        ('07_pairrule_peak_screen', 'figure7_ver05_analysis/pairrule', [py, 'build_peak_screen.py']),
        ('08_gap_peak_extraction', 'figure7_ver05_analysis/gap_peaks', [py, 'extract_gap_peaks.py']),
        ('09_all_gene_comparison', 'figure7_ver05_analysis/gap_peaks', [py, 'compare_peak_orders.py']),
        ('10_s1_data', 'figure7_ver05', [py, 'assemble_outputs.py']),
    ]
    if not args.skip_figures:
        steps.extend([
            ('11_figure7', 'figure7_ver08', [py, 'build_figure7.py']),
            ('12_companion', 'figure7_ver06', [py, 'build_companion.py']),
        ])
    timings = []; started = time.monotonic()
    for label, folder, command in steps:
        print(f'Running {label}...', flush=True)
        tick = time.monotonic()
        with (logdir / f'{label}.txt').open('w') as log:
            run = subprocess.run(command, cwd=work / folder, stdout=log, stderr=subprocess.STDOUT)
        if run.returncode:
            raise RuntimeError(f'{label} failed. See {logdir / (label + ".txt")}')
        timings.append({'step': label, 'seconds': round(time.monotonic()-tick, 3)})
    # Named submission outputs, without requiring readers to know historical paths.
    shutil.copyfile(work / 'figure7_ver05/Peak_Order_Comparisons_ver05.csv', output / 'S1_Data_Peak_Order_Comparisons.csv')
    shutil.copyfile(work / 'figure7_ver05/Peak_Order_Gene_Pair_Coverage_ver05.csv', output / 'Gene_Pair_Coverage.csv')
    shutil.copyfile(work / 'figure7_ver05_analysis/pairrule/posterior_mask.csv', output / 'Posterior_Domain_Mask.csv')
    if not args.skip_figures:
        for ext in ['pdf', 'svg', 'png', 'json']:
            shutil.copyfile(work / f'figure7_ver08/Fig7_ver08.{ext}', output / f'Figure_7.{ext}')
        shutil.copyfile(work / 'figure7_ver06/Fig7_Peak_Order_Companion_ver06.pdf', output / 'Peak_Order_Companion.pdf')
        for name in ['Fig7_peak_time_differences_ver08.csv', 'Fig7_peak_markers_ver08.csv']:
            shutil.copyfile(work / 'figure7_ver08' / name, output / name)
    validation = []
    if args.bootstraps == 20000 and not args.no_reference_check and not args.skip_figures and not args.apply_confirmed_corrections:
        validation = validate_tables(work)
    summary = {'bootstrap_draws': args.bootstraps, 'seed': 20261019,
               'confirmed_temporal_corrections_applied': args.apply_confirmed_corrections,
               'seconds': round(time.monotonic()-started, 3), 'steps': timings,
               'reference_checks': validation, 'python': sys.version,
               'scientific_input_sha256': {str(p.relative_to(PACKAGE)): hashlib.sha256(p.read_bytes()).hexdigest()
                                           for p in sorted((PACKAGE / 'inputs').rglob('*')) if p.is_file()}}
    summary['repository_input_sha256'] = {str(p.relative_to(destination)): hashlib.sha256(p.read_bytes()).hexdigest()
                                        for p in sorted(destination.rglob('*')) if p.is_file()}
    (output / 'run_manifest.json').write_text(json.dumps(summary, indent=2))
    print(f'Done. Results: {output}', flush=True)
    print(f'Reference tables passed: {len(validation)}', flush=True)


if __name__ == '__main__':
    main()
