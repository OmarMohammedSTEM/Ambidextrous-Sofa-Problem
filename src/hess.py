# Full finite-difference Hessian of A(cR, cL) at the symmetric K=4 optimum,
# in the 36-dim asymmetric (both-paths-free) space. Central differences.
import json, time, numpy as np
import ambi

r = json.load(open('data/ambi_results.json'))
x = np.array(r['B_K4']['coef'])          # 18 params
n = len(x); m = 2*n                      # 36-dim (cR, cL)
z0 = np.concatenate([x, x])
N = 161; h = 0.015

def f(z):
    return ambi.area(z[:n], z[n:], N)

t0 = time.time()
base = f(z0)
H = np.zeros((m, m))
# diagonal
fp = np.zeros(m); fm = np.zeros(m)
for i in range(m):
    zp = z0.copy(); zp[i] += h
    zm = z0.copy(); zm[i] -= h
    fp[i] = f(zp); fm[i] = f(zm)
    H[i, i] = (fp[i] + fm[i] - 2*base)/h**2
print(f'diagonal done [{time.time()-t0:.0f}s]', flush=True)
# off-diagonal
for i in range(m):
    for j in range(i+1, m):
        zpp = z0.copy(); zpp[i] += h; zpp[j] += h
        zpm = z0.copy(); zpm[i] += h; zpm[j] -= h
        zmp = z0.copy(); zmp[i] -= h; zmp[j] += h
        zmm = z0.copy(); zmm[i] -= h; zmm[j] -= h
        H[i, j] = H[j, i] = (f(zpp) - f(zpm) - f(zmp) + f(zmm))/(4*h**2)
    if i % 6 == 0:
        print(f'row {i}/{m} [{time.time()-t0:.0f}s]', flush=True)

evals, evecs = np.linalg.eigh(H)
# gauge direction: joint horizontal translation = a0 of both paths together
g = np.zeros(m); g[0] = g[n] = 1/np.sqrt(2)
ovl = evecs.T @ g
top = np.argmax(evals)
print('eigenvalues (ascending):', np.array2string(evals, precision=3), flush=True)
print(f'largest eigenvalue {evals[top]:+.4f}; overlap of its eigenvector with translation gauge: {abs(ovl[top]):.3f}', flush=True)
print(f'#eig < -0.1: {(evals < -0.1).sum()}, in [-0.1, -0.005]: {((evals >= -0.1)&(evals < -0.005)).sum()}, in [-0.005, 0.005]: {(abs(evals) <= 0.005).sum()}, > 0.005: {(evals > 0.005).sum()}', flush=True)
json.dump(dict(eigenvalues=list(evals), gauge_overlap=list(map(float, abs(ovl))),
               h=h, N=N, base=base), open('data/hess.json', 'w'), indent=1)
print(f'DONE [{time.time()-t0:.0f}s]', flush=True)
