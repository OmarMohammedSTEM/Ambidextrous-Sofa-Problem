# Hessian at parametrized step size h (usage: hess2.py <h>), saving H and eigenvectors.
# Purpose: h-stability check of the soft eigenvalues vs the gauge-calibrated noise floor.
import json, sys, time, numpy as np
import ambi

h = float(sys.argv[1])
r = json.load(open('data/ambi_results.json'))
x = np.array(r['B_K4']['coef'])
n = len(x); m = 2*n
z0 = np.concatenate([x, x])
N = 161

def f(z): return ambi.area(z[:n], z[n:], N)

t0 = time.time()
base = f(z0)
H = np.zeros((m, m))
for i in range(m):
    zp = z0.copy(); zp[i] += h; zm = z0.copy(); zm[i] -= h
    H[i, i] = (f(zp) + f(zm) - 2*base)/h**2
for i in range(m):
    for j in range(i+1, m):
        zpp = z0.copy(); zpp[i] += h; zpp[j] += h
        zpm = z0.copy(); zpm[i] += h; zpm[j] -= h
        zmp = z0.copy(); zmp[i] -= h; zmp[j] += h
        zmm = z0.copy(); zmm[i] -= h; zmm[j] -= h
        H[i, j] = H[j, i] = (f(zpp) - f(zpm) - f(zmp) + f(zmm))/(4*h**2)
    if i % 9 == 0: print(f'h={h}: row {i}/{m} [{time.time()-t0:.0f}s]', flush=True)

ev, V = np.linalg.eigh(H)
g = np.zeros(m); g[0] = g[n] = 1/np.sqrt(2)
ovl = abs(V.T @ g)
print(f'h={h}: top 5 eigenvalues {np.round(ev[-5:], 4)}; gauge overlap of top-3: {np.round(ovl[-3:],3)}', flush=True)
json.dump(dict(h=h, N=N, eigenvalues=list(ev), gauge_overlap=list(map(float, ovl)),
               H=[list(row) for row in H]),
          open(f'data/hess_h{h}.json', 'w'))
print(f'DONE [{time.time()-t0:.0f}s]', flush=True)
