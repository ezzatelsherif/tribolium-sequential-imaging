"""Regenerate the final-fit mesh refinement table without nonlinear refitting."""
import json
import numpy as np
import pandas as pd
from pulse_fit import ROOT, GENES, PulseFit


def main():
    rows = []
    for gene in GENES:
        for shape in ('cosine', 'gaussian'):
            fits = json.loads((ROOT / 'fitresults' / f'fit_{gene}_{shape}.json').read_text())
            for mode, fit in fits.items():
                y, original = np.asarray(fit['y']), np.asarray(fit['pred'])
                for dt in (.025, .0125, .00625):
                    model = PulseFit(gene, mode=mode, shape=shape, dt=dt)
                    prediction = model.predict(fit)
                    profiled = model.solve(model.pack(fit['parameters']), y)
                    rows.append(dict(gene=gene, shape=shape, mode=mode, dt=dt,
                        SSE_fixed_coefficients=float(np.sum((prediction-y)**2)),
                        SSE_reprofiled_coefficients=profiled['cost'],
                        maximum_prediction_change=float(np.max(np.abs(prediction-original)))))
    output = pd.DataFrame(rows)
    output.to_csv(ROOT / 'fitresults' / 'mesh_refinement.csv', index=False)
    print(output.groupby('dt').maximum_prediction_change.max().to_string())


if __name__ == '__main__':
    main()
