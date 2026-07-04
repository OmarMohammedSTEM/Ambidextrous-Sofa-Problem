# K=5 continuation: does the deficit to Romik keep shrinking with bandwidth?
import json, time, numpy as np
import ambi
from scipy.optimize import minimize

r = json.load(open('data/ambi_results.json'))
x = ambi.pad(np.array(r['B_K4']['coef']), 5)   # 22 params
t0 = time.time()
res = minimize(lambda p: -ambi.area(p, p, 101), x, method='Powell',
               options=dict(maxfev=4500, xtol=1e-7, ftol=1e-10))
x = res.x
ext, a1, a2 = ambi.fine_area(x, x, N=321)
print(f'K=5: raw101={-res.fun:.7f}  A(321)={a1:.7f} A(641)={a2:.7f} Richardson={ext:.7f} '
      f'deficit={1.644955218-ext:+.3e}  nfev={res.nfev} [{time.time()-t0:.0f}s]', flush=True)
r['B_K5'] = dict(coef=list(x), area=ext, nfev=int(res.nfev))
json.dump(r, open('data/ambi_results.json', 'w'), indent=1)
print('DONE', flush=True)
