"""Portable entry points for the intron/exon analysis.

The default verifies saved results without optimization. The all stage recomputes
empirical estimates, the hypothetical sweep, fits, assessment, and numerical checks.
"""
from pathlib import Path
import argparse
import json
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def run(*args):
    subprocess.run([sys.executable, *map(str, args)], cwd=ROOT, check=True)


def main(stage, fresh=False):
    if stage == 'verify':
        run('verify_saved_results.py')
        return
    if stage in ('manuscript-figures', 'movie'):
        extra = ['--movie'] if stage == 'movie' else []
        run('make_manuscript_figures.py', '--numerical-only', *extra)
        return
    if stage in ('all', 'empirical'):
        derived = ROOT / 'empirical' / 'derived'
        run('empirical/original_stage_lag_analysis.py', '--output-dir', derived)
        for name in ('stage_lag_bootstrap.csv', 'stage_lag_sensitivity.csv'):
            shutil.copy2(derived / name, ROOT / 'empirical' / name)
        shutil.copy2(derived / 'background_corrected_pairs.csv',
                     ROOT / 'data' / 'background_corrected_pairs.csv')
    if stage in ('all', 'sweep'):
        run('sweep.py')
    if stage in ('all', 'fits'):
        for shape in ('cosine', 'gaussian'):
            for gene in ('Kr', 'mlpt', 'svb'):
                extra = ['--fresh'] if fresh else []
                run('pulse_fit.py', '--gene', gene, '--shape', shape, *extra)
    if stage in ('all', 'numerics'):
        run('check_fit_mesh.py')
        run('validate_pulse_numerics.py')
    if stage in ('all', 'assessment'):
        for gene in ('Kr', 'mlpt', 'svb'):
            run('assess_constrained.py', '--gene', gene, '--task', 'profile')
            run('assess_constrained.py', '--gene', gene, '--task', 'null', '--draws', '199')
            run('assess_constrained.py', '--gene', gene, '--task', 'check')
        run('assess_constrained.py', '--task', 'summary')
    if stage in ('all', 'figures'):
        run('make_sweep_figures.py')
        run('make_fit_figures.py')
        run('make_empirical_figure.py')
    if stage == 'all':
        run('verify_saved_results.py', '--skip-manifest')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['verify', 'all', 'empirical', 'sweep',
        'fits', 'numerics', 'assessment', 'figures', 'manuscript-figures', 'movie'], default='verify')
    parser.add_argument('--fresh', action='store_true',
        help='For fits/all, ignore archived fit parameters as warm starts.')
    args = parser.parse_args()
    main(args.stage, args.fresh)
