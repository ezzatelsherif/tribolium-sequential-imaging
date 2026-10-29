#!/usr/bin/env python3
"""Rebuild the landmark inventory and pooled summary using the saved biological draws."""
from pathlib import Path
import subprocess,sys,json,hashlib,argparse,shutil
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figure',action='store_true',help='Also rebuild the final Figure 7 artwork; requires requirements-figure.txt and Inkscape on PATH.')
    args=parser.parse_args()
    subprocess.run([sys.executable,str(ROOT/'build_dataset.py')],cwd=ROOT,check=True)
    subprocess.run([sys.executable,str(ROOT/'summarize_landmarks.py')],cwd=ROOT,check=True)
    comparisons={}
    for generated,reference in [
        ('data/landmark_distances.csv','reference/landmark_distances.csv'),
        ('outputs/S2_Data_Landmark_Comparisons.csv','reference/S2_Data_Landmark_Comparisons.csv'),
        ('outputs/S3_Data_Landmark_Summary.csv','reference/S3_Data_Landmark_Summary.csv'),
    ]:
        a,b=ROOT/generated,ROOT/reference
        if not b.exists():
            raise FileNotFoundError(b)
        exact=a.read_bytes()==b.read_bytes()
        if not exact:
            pd.testing.assert_frame_equal(pd.read_csv(a),pd.read_csv(b),check_dtype=False,rtol=0,atol=1e-12)
        comparisons[generated]={'byte_identical':exact,'numerical_match':True,'sha256':sha(a)}
    features=np.load(ROOT/'data/temporal_features_with_observed.npz')
    contrasts=np.load(ROOT/'data/temporal_contrast_draws.npz')
    assert all(features[k].shape==(20001,) for k in features.files)
    assert all(contrasts[k].shape==(20000,) for k in contrasts.files)
    report={'reconstruction':'saved biological draws; no embryo resampling rerun','comparisons':comparisons,'feature_arrays':len(features.files),'contrast_arrays':len(contrasts.files),'bootstrap_draws':20000,'all_checks_pass':True}
    workbook = ROOT/'reference/S1_Data_Landmark_Results.xlsx'
    subprocess.run([sys.executable,str(ROOT/'verify_results_workbook.py'),
                    '--workbook',str(workbook),
                    '--peak-audit',str(ROOT/'data/source_tables/peak_inventory.csv'),
                    '--landmark-inventory',str(ROOT/'outputs/S2_Data_Landmark_Comparisons.csv'),
                    '--summary-source',str(ROOT/'outputs/S3_Data_Landmark_Summary.csv'),
                    '--report',str(ROOT/'validation/workbook_replay.json')],cwd=ROOT,check=True)
    shutil.copyfile(workbook,ROOT/'outputs/S1_Data_Landmark_Results.xlsx')
    report['workbook_verified_against_reconstructed_tables']=True
    report['workbook_sha256']=sha(workbook)
    if args.figure:
        figure=ROOT/'figure7'
        shutil.copyfile(ROOT/'outputs/S3_Data_Landmark_Summary.csv',figure/'source/summary_counts_and_percentages.csv')
        shutil.copyfile(ROOT/'outputs/S2_Data_Landmark_Comparisons.csv',figure/'source/S2_Data_Landmark_Comparisons.csv')
        subprocess.run([sys.executable,str(figure/'build_fig7_ver12.py')],cwd=ROOT,check=True)
        assert (figure/'Fig7_ver12.svg').read_bytes()==(figure/'reference/Fig7_ver12.svg').read_bytes(), 'Rebuilt figure SVG differs from the frozen reference'
        report['figure_svg_byte_identical']=True
        report['figure_validation']=json.loads((figure/'qa/validation.json').read_text())
    (ROOT/'validation').mkdir(exist_ok=True)
    (ROOT/'validation/replay.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
