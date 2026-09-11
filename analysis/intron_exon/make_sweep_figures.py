"""Replot the recorded sweep without rerunning model scenarios."""
import pandas as pd
from sweep import ROOT, make_figures

if __name__ == '__main__':
    make_figures(pd.read_csv(ROOT / 'sweep_results' / 'sweep_results.csv'))
