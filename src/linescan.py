# Direct line scans of the area along the softest Hessian eigenvectors (Section 7).
# Bypasses finite differences: gauge direction must return machine zero (exact
# translation invariance); genuine soft modes must decrease ~ (1/2) lambda a^2.
# Requires hess_h0.03.json (full H matrix) from `hess2.py 0.03`.
import json, numpy as np
import ambi

d = json.load(open('data/hess_h0.03.json'))
H = np.array(d['H'])
ev, V = np.linalg.eigh(H)
x = np.array(json.load(open('data/ambi_results.json'))['B_K4']['coef']); n = len(x)
z0 = np.concatenate([x, x])
N = 161
base = ambi.area(x, x, N)
g = np.zeros(2*n); g[0] = g[n] = 1/np.sqrt(2)

print('top6 eigenvalues (h=0.03):', np.round(ev[-6:], 4))
results = {}
for k in range(1, 7):
    v = V[:, -k]
    row = []
    for amp in (0.01, 0.02, 0.05, 0.1):
        z = z0 + amp*v
        row.append(ambi.area(z[:n], z[n:], N) - base)
    results[f'mode_{-k}'] = dict(lam=float(ev[-k]), gauge_ovl=float(abs(v @ g)),
                                 amps=[0.01, 0.02, 0.05, 0.1], dA=row)
    print(f'mode {-k} (lam={ev[-k]:+.4f}, gauge_ovl={abs(v @ g):.2f}): '
          + '  '.join(f'dA({a})={dA:+.1e}' for a, dA in zip((0.01,0.02,0.05,0.1), row)))
    # quadratic fit dA = (1/2) lam_eff a^2 from the two largest amplitudes
    lam_eff = 2*row[-1]/0.1**2
    print(f'         effective lambda from a=0.1: {lam_eff:+.4f}')

json.dump(results, open('data/linescan.json', 'w'), indent=1)
print('DONE -> linescan.json')
