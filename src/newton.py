# Newton polish at the K=5 point: if the deficit is optimizer non-convergence,
# a curvature-corrected step should reduce it measurably.
import json, time, numpy as np
import ambi

r = json.load(open('data/ambi_results.json'))
x = np.array(r['B_K5']['coef']); n = len(x)
N = 161; t0 = time.time()

def A(p): return ambi.area(p, p, N)

for rnd in range(2):
    base = A(x)
    hg = 0.005
    g = np.zeros(n)
    for i in range(n):
        xp = x.copy(); xm = x.copy(); xp[i] += hg; xm[i] -= hg
        g[i] = (A(xp) - A(xm))/(2*hg)
    hh = 0.015
    H = np.zeros((n, n))
    for i in range(n):
        xp = x.copy(); xm = x.copy(); xp[i] += hh; xm[i] -= hh
        H[i, i] = (A(xp) + A(xm) - 2*base)/hh**2
    for i in range(n):
        for j in range(i+1, n):
            zpp = x.copy(); zpp[i] += hh; zpp[j] += hh
            zpm = x.copy(); zpm[i] += hh; zpm[j] -= hh
            zmp = x.copy(); zmp[i] -= hh; zmp[j] += hh
            zmm = x.copy(); zmm[i] -= hh; zmm[j] -= hh
            H[i, j] = H[j, i] = (A(zpp) - A(zpm) - A(zmp) + A(zmm))/(4*hh**2)
    ev, V = np.linalg.eigh(H)
    print(f'round {rnd}: base={base:.7f} |g|={np.linalg.norm(g):.3e} '
          f'eig range [{ev[0]:.2f}, {ev[-1]:.3f}] [{time.time()-t0:.0f}s]', flush=True)
    # Newton direction restricted to well-conditioned negative subspace
    d = np.zeros(n)
    for lam, v in zip(ev, V.T):
        if lam < -0.05:
            d -= (v @ g)/lam * v
    if np.linalg.norm(d) > 0.08: d *= 0.08/np.linalg.norm(d)
    best_a, best_alpha = base, 0.0
    for alpha in (0.25, 0.5, 1.0, 1.5):
        a = A(x + alpha*d)
        print(f'  alpha={alpha}: A={a:.7f} ({a-base:+.2e})', flush=True)
        if a > best_a: best_a, best_alpha = a, alpha
    if best_alpha == 0.0:
        print('no improvement; stopping', flush=True); break
    x = x + best_alpha*d

ext, a1, a2 = ambi.fine_area(x, x, N=321)
print(f'after Newton: A(321)={a1:.7f} A(641)={a2:.7f} Richardson={ext:.7f} '
      f'deficit={1.644955218-ext:+.3e} [{time.time()-t0:.0f}s]', flush=True)
r['B_K5_newton'] = dict(coef=list(x), area=ext)
json.dump(r, open('data/ambi_results.json', 'w'), indent=1)
print('DONE', flush=True)
