"""
Free-knot cubic spline pathspace — kills the algebraic Fourier tail.

Romik's exact corner path has curvature jumps at contact-switch angles
(approx 0.3 and 1.27 rad for the ambidextrous case; Gerver's single-turn has
well-known kinks at similar positions). A global Fourier basis converges
only as O(1/k^2) there, which is why K=5 still leaves ~1e-5 deficit.

A cubic spline with knots *at* those angles represents the kinks exactly
with O(h^4) smooth-piece accuracy. With free knots we let the optimizer
place them.

Parametrization:
  xi(t) = CubicSpline(t_knots, x_nodes)
  eta(t) = CubicSpline(t_knots, y_nodes)
  knots = [0, k1, k2, pi/2 - k2, pi/2 - k1, pi/2]  (symmetric about pi/4)
  k1, k2 in (0, pi/4) are free variables, ordered 0 < k1 < k2 < pi/4.
"""

import json, time, sys, pathlib
import numpy as np

# Ensure anc importable with -I
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

from scipy.interpolate import CubicSpline
from scipy.optimize import minimize

PI2 = np.pi/2
ROMIK = 1.64495521842544

def build_knots(k1, k2):
    k1 = float(np.clip(k1, 0.05, 0.70))
    k2 = float(np.clip(k2, 0.10, 0.75))
    if k1 >= k2: k1, k2 = k2, k1+1e-6
    if k2 >= PI2/2: k2 = PI2/2 - 0.02
    knots = np.array([0.0, k1, k2, PI2 - k2, PI2 - k1, PI2])
    knots = np.sort(knots)
    return knots

def eval_spline_area(x_nodes, y_nodes, knots, N):
    cs_x = CubicSpline(knots, x_nodes, bc_type='natural')
    cs_y = CubicSpline(knots, y_nodes, bc_type='natural')
    ts = np.linspace(0, PI2, N)
    xi = cs_x(ts)
    eta = cs_y(ts)
    from shapely.geometry import Polygon
    from shapely import affinity
    S = ambi.STRIP
    H0 = ambi.H0
    for t, a, b in zip(ts, xi, eta):
        Ht = affinity.translate(affinity.rotate(H0, float(t), origin=(0,0), use_radians=True), float(a), float(b))
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

def optimize_spline(N_search=81, N_fine=321, iters=2500, seed_knots=(0.35, 0.85)):
    k1_0, k2_0 = seed_knots
    knots0 = build_knots(k1_0, k2_0)
    m = len(knots0)

    # Warm start: sample the B_K5_newton Fourier path onto the spline nodes
    xi0 = eta0 = None
    for _p in ('anc/ambi_results.json', 'ambi_results.json', 'C:/New folder (31)/anc/ambi_results.json'):
        try:
            if not pathlib.Path(_p).exists():
                continue
            d = json.load(open(_p))
            for _k in ('B_K5_newton', 'B_K5', 'B_K4'):
                if _k in d and 'coef' in d[_k]:
                    coef = np.array(d[_k]['coef'])
                    xi0, eta0 = ambi.make_path(coef, knots0)
                    break
            if xi0 is not None:
                break
        except Exception:
            continue
    if xi0 is None:
        r = 0.325
        xi0 = r*np.cos(2*knots0)
        eta0 = r*np.sin(2*knots0)
        print('[spline] warm start fallback: Hammersley circle', flush=True)

    x0 = np.concatenate([xi0, eta0, np.array([k1_0, k2_0])])

    def neg_area(p):
        m2 = m
        x_nodes = p[0:m2]
        y_nodes = p[m2:2*m2]
        k1, k2 = float(p[-2]), float(p[-1])
        k1 = np.clip(k1, 0.08, 0.70)
        k2 = np.clip(k2, 0.15, 0.75)
        if k1 >= k2:
            return 1e6
        knots = build_knots(k1, k2)
        return -eval_spline_area(x_nodes, y_nodes, knots, N_search)

    print(f'[spline] m={m} knots, N_search={N_search}, Powell maxfev={iters}', flush=True)
    print(f'[spline] initial knots {np.round(knots0,5)}', flush=True)
    a0 = -neg_area(x0)
    print(f'[spline] initial area N={N_search}: {a0:.6f}', flush=True)
    t0 = time.time()
    res = minimize(neg_area, x0, method='Powell',
                   options=dict(maxfev=iters, xtol=1e-7, ftol=1e-10))
    p_opt = res.x
    m2 = m
    x_opt = p_opt[0:m2]; y_opt = p_opt[m2:2*m2]
    k1_opt, k2_opt = float(p_opt[-2]), float(p_opt[-1])
    knots_opt = build_knots(k1_opt, k2_opt)
    print(f'[spline] Powell done nfev={res.nfev} [{time.time()-t0:.0f}s]', flush=True)
    print(f'[spline] optimal knots {np.round(knots_opt,5)}', flush=True)

    def fine(knots, xn, yn):
        a1 = eval_spline_area(xn, yn, knots, N_fine)
        a2 = eval_spline_area(xn, yn, knots, 2*N_fine-1)
        return 2*a2 - a1, a1, a2

    ext, a1, a2 = fine(knots_opt, x_opt, y_opt)
    print(f'[spline] fine N={N_fine}: A={a1:.7f} A2={a2:.7f} Richardson={ext:.7f}', flush=True)
    print(f'[spline] deficit to Romik: {ROMIK - ext:+.3e}', flush=True)

    out = dict(knots=list(map(float, knots_opt)), x_nodes=list(map(float, x_opt)),
               y_nodes=list(map(float, y_opt)), k1=float(k1_opt), k2=float(k2_opt),
               area_richardson=float(ext), area_N1=float(a1), area_N2=float(a2),
               nfev=int(res.nfev))
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
    prev['spline_free_knot'] = out
    json.dump(prev, open(save_path,'w'), indent=1)
    print(f'[spline] saved -> {save_path}["spline_free_knot"]', flush=True)
    return out

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=81)
    ap.add_argument('--iters', type=int, default=2500)
    args = ap.parse_args()
    optimize_spline(N_search=args.N, iters=args.iters)
