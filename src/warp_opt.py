"""
Non-monotone rotation test.

Section 9 caveat: the Fourier basis fixes rotation schedule theta(t)=t
(monotone). Could a non-monotone angle schedule — revisiting an angle or
overshooting pi/2 within a single turn — improve the ambidextrous area?

Heuristic: revisiting an angle adds a *superset* of constraints
(intersection over T including repeated angles is subset of monotone intersection),
so non-monotone should be dominated. But the corridor stations for the two
turns are distinct, so the heuristic is not a proof; the second family's
mirror could in principle compensate.

We test numerically by parametrizing the rotation schedule itself:

  theta(t) = t + sum_{k=1}^{Km}  e_k sin(2 k t)    (deviation from identity)
  with theta(0)=0, theta(pi/2)=pi/2 enforced by the sine basis (zero at endpoints).

If all e_k = 0 we recover monotone. Non-monotone corresponds to regions where
theta'(t) = 1 + 2k e_k cos(2kt) < 0 for some t.

We jointly optimize (Fourier path coefficients c) and (angle-warp coefficients e)
to maximize area. If the optimum stays at e approx 0 (or area does not improve), that
is strong evidence that monotone is optimal within smooth schedules.

We use the same Shapely kernel: instead of uniform t_i, we use theta_i = theta(t_i)
as the hallway rotation angles while the corner still rides v(t_i).
"""

import json, time, sys, pathlib
import numpy as np
from scipy.optimize import minimize

# Ensure anc is importable even with python -I
for _cand in ('anc', '.', 'C:/New folder (31)/anc', 'C:/New folder (31)'):
    if _cand not in sys.path and pathlib.Path(_cand).exists():
        sys.path.insert(0, _cand)
try:
    import ambi
except Exception:
    import importlib.util as _ilu
    for _p in ('anc/ambi.py', 'C:/New folder (31)/anc/ambi.py', 'ambi.py'):
        if pathlib.Path(_p).exists():
            _spec = _ilu.spec_from_file_location('ambi', _p)
            ambi = _ilu.module_from_spec(_spec); _spec.loader.exec_module(ambi)
            break

PI2 = np.pi/2
ROMIK = 1.64495521842544

def theta_warp(t, e):
    t = np.asarray(t, float)
    th = t.copy()
    for k, ek in enumerate(e, start=1):
        th = th + ek * np.sin(2*k*t)
    return th

def theta_of_t(t, e_coef):
    return theta_warp(t, e_coef)

def theta_warp_deriv(t, e):
    t = np.asarray(t, float)
    d = np.ones_like(t)
    for k, ek in enumerate(e, start=1):
        d = d + ek * 2*k * np.cos(2*k*t)
    return d

def area_with_warp(coef, e, N):
    """Sofa area when hallway i is rotated by theta_i = theta_warp(t_i, e)
    and translated to v(t_i). Both families use the same warp (symmetric case).
    For a fully general test we allow independent warps eR,eL — here we keep
    one warp for the symmetric section."""
    ts = np.linspace(0, PI2, N)
    xi, eta = ambi.make_path(coef, ts)
    th = theta_warp(ts, e)
    from shapely.geometry import Polygon
    from shapely import affinity
    S = ambi.STRIP
    H0 = ambi.H0
    for t_raw, th_i, a, b in zip(ts, th, xi, eta):
        Ht = affinity.translate(affinity.rotate(H0, float(th_i), origin=(0,0), use_radians=True), float(a), float(b))
        S = S.intersection(Ht)
        if S.is_empty:
            return 0.0
        Hm = affinity.scale(Ht, 1.0, -1.0, origin=(0.0, 0.5))
        S = S.intersection(Hm)
        if S.is_empty:
            return 0.0
    if S.geom_type != 'Polygon':
        polys = [g for g in getattr(S, 'geoms', [S]) if g.geom_type == 'Polygon']
        if not polys:
            return 0.0
        S = max(polys, key=lambda p: p.area)
    return S.area

def optimize_warp(Km=2, N_search=61, maxfev=2500):
    # Load seed coef robustly
    coef0 = None
    for _p in ('anc/ambi_results.json', 'ambi_results.json', 'C:/New folder (31)/anc/ambi_results.json'):
        try:
            if not pathlib.Path(_p).exists():
                continue
            d = json.load(open(_p))
            for _k in ('B_K5_newton', 'B_K5', 'B_K4', 'B_K3'):
                if _k in d and 'coef' in d[_k]:
                    coef0 = np.array(d[_k]['coef']); break
            if coef0 is not None:
                break
        except Exception:
            continue
    if coef0 is None:
        coef0 = np.array([0,0,0.325,0,0,0.325], float)

    e0 = np.zeros(Km)
    z0 = np.concatenate([coef0, e0])
    n_coef = len(coef0)

    def neg_area(z):
        coef = z[:n_coef]
        e = z[n_coef:n_coef+Km]
        th = theta_warp(np.linspace(0, PI2, 9), e)
        if np.any(th < -0.25) or np.any(th > PI2+0.25):
            return 1e6
        return -area_with_warp(coef, e, N_search)

    print(f'[warp] Km={Km} N={N_search} joint Powell over {len(z0)} vars (coef {n_coef} + warp {Km})', flush=True)
    base = -neg_area(z0)
    print(f'[warp] baseline (e=0) area N={N_search}: {base:.6f}', flush=True)
    t0 = time.time()
    res = minimize(neg_area, z0, method='Powell',
                   options=dict(maxfev=maxfev, xtol=1e-7, ftol=1e-10))
    coef_opt = res.x[:n_coef]
    e_opt = res.x[n_coef:n_coef+Km]
    best = -res.fun
    print(f'[warp] Powell done nfev={res.nfev} [{time.time()-t0:.0f}s]', flush=True)
    print(f'[warp] e_opt = {np.round(e_opt,5)}  max|th - t| = {np.max(np.abs(theta_warp(np.linspace(0,PI2,200), e_opt)-np.linspace(0,PI2,200))):.4f} rad', flush=True)
    ts_dense = np.linspace(0, PI2, 400)
    deriv = theta_warp_deriv(ts_dense, e_opt)
    print(f'[warp] min theta_prime(t) = {deriv.min():+.4f}  (negative => non-monotone)', flush=True)
    print(f'[warp] best area N={N_search}: {best:.6f}  delta vs baseline {best-base:+.2e}', flush=True)
    def fine_warp(coef, e):
        a1 = area_with_warp(coef, e, 241)
        a2 = area_with_warp(coef, e, 481)
        return 2*a2 - a1, a1, a2
    ext, a1, a2 = fine_warp(coef_opt, e_opt)
    ext0, _, _ = fine_warp(coef0, np.zeros(Km))
    print(f'[warp] Richardson fine: warped {ext:.7f} vs unwarped {ext0:.7f}  (Romik {ROMIK:.7f})', flush=True)

    out = dict(Km=Km, e_opt=list(map(float, e_opt)), coef_opt=list(map(float, coef_opt)),
               area_warped=float(ext), area_unwarped=float(ext0), best_N_area=float(best),
               min_deriv=float(deriv.min()), nfev=int(res.nfev))
    # save
    save_path = None
    for _p in ('anc/ambi_results.json', 'C:/New folder (31)/anc/ambi_results.json'):
        if pathlib.Path(_p).exists():
            save_path = _p; break
    if save_path is None:
        save_path = 'anc/ambi_results.json'
    try:
        prev = json.load(open(save_path))
    except Exception:
        prev = {}
    prev['warp_test'] = out
    json.dump(prev, open(save_path,'w'), indent=1)
    print(f'[warp] saved -> {save_path}["warp_test"]', flush=True)
    return out

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--Km', type=int, default=2)
    ap.add_argument('--N', type=int, default=61)
    ap.add_argument('--iters', type=int, default=2500)
    args = ap.parse_args()
    optimize_warp(Km=args.Km, N_search=args.N, maxfev=args.iters)
