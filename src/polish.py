import json, time, numpy as np
import ambi
from scipy.optimize import minimize

r = json.load(open('data/ambi_results.json'))
x = ambi.pad(np.array(r['B_K3']['coef']), 4)

t0 = time.time()
f = lambda p: -ambi.area(p, p, 101)
res = minimize(f, x, method='Powell', options=dict(maxfev=3500, xtol=1e-7, ftol=1e-10))
x = res.x
ext, a1, a2 = ambi.fine_area(x, x, N=321)
print(f'K=4 polish: raw101={-res.fun:.7f}  Richardson(321/641)={ext:.7f} '
      f'deficit vs Romik {1.644955218-ext:+.2e}  nfev={res.nfev} [{time.time()-t0:.0f}s]', flush=True)
r['B_K4'] = dict(coef=list(x), area=ext)
json.dump(r, open('data/ambi_results.json','w'), indent=1)

# ---- figure ----
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
S = ambi.sofa_poly(x, x, 641)
xs, ys = np.array(S.exterior.xy)
xs -= (xs.max()+xs.min())/2  # center
ts = np.linspace(0, np.pi/2, 400)
xi, eta = ambi.make_path(x, ts)
xi -= np.mean([xi.max()+xi.min()])/1e9  # keep as-is; center via same shift as sofa
# recompute proper shift: use uncentered polygon extent
S2 = ambi.sofa_poly(x, x, 641)
ux, uy = np.array(S2.exterior.xy)
cx = (ux.max()+ux.min())/2

fig, ax = plt.subplots(figsize=(9,4))
ax.fill(xs, ys, color='#cfe3f7', ec='#1f5fa8', lw=1.5, zorder=2,
        label=f'optimized ambidextrous sofa  A={ext:.5f}')
ax.plot(xi-cx, eta, 'r--', lw=1.2, label='corner path v(t): lower notch (right turn)')
ax.plot(xi-cx, 1-eta, 'g--', lw=1.2, label='mirror path (left turn)')
ax.axhline(0, color='k', lw=.6); ax.axhline(1, color='k', lw=.6)
ax.set_aspect('equal'); ax.legend(loc='center', fontsize=8, framealpha=.9)
ax.set_title("Emergent ambidextrous sofa (Fourier path space, K=4)  vs  Romik's 1.6449552")
plt.tight_layout(); plt.savefig('figures/ambi_sofa.png', dpi=150)
print('figure saved: ambi_sofa.png', flush=True)
print('final coefficients:', ['%.5f'%c for c in x], flush=True)
