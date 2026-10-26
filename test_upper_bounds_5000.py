import sys, time, json, numpy as np
from shapely.geometry import Polygon
from shapely import affinity
from scipy.interpolate import interp1d

ARM = 6.0
H0 = Polygon([(1,1),(-ARM,1),(-ARM,0),(0,0),(0,-ARM),(1,-ARM)])
STRIP = Polygon([(-ARM,0),(ARM,0),(ARM,1),(-ARM,1)])

def discrete_sofa_area(vR, wL, tR, tL):
    S = STRIP
    for t, p in zip(tR, vR):
        Ht = affinity.translate(
            affinity.rotate(H0, float(t), origin=(0,0), use_radians=True),
            float(p[0]), float(p[1]))
        S = S.intersection(Ht)
        if S.is_empty: return 0.0
    for s, q in zip(tL, wL):
        Hs = affinity.translate(
            affinity.rotate(H0, float(s), origin=(0,0), use_radians=True),
            float(q[0]), float(q[1]))
        Hs = affinity.scale(Hs, 1.0, -1.0, origin=(0.0, 0.5))
        S = S.intersection(Hs)
        if S.is_empty: return 0.0
    if S.geom_type != 'Polygon':
        polys = [g for g in getattr(S, 'geoms', [S]) if g.geom_type == 'Polygon']
        if not polys: return 0.0
        S = max(polys, key=lambda p: p.area)
    return S.area

# Load N=30 solution as base
d = json.load(open('C:/New folder (31)/upper_bound_results.json'))
t_base = np.array(d['30']['tR'])
vR_base = np.array(d['30']['vR'])
wL_base = np.array(d['30']['wL'])

ROMIK = 1.644955218

print("=== Pushing Upper Bounds to N=5000 (20,000 variables) ===", flush=True)
N = 5000
t0 = time.time()
t_new = np.linspace(0, np.pi/2, N)
f_vR = interp1d(t_base, vR_base, axis=0, kind='cubic', fill_value='extrapolate')
f_wL = interp1d(t_base, wL_base, axis=0, kind='cubic', fill_value='extrapolate')
vR_init = f_vR(t_new)
wL_init = f_wL(t_new)

a_init = discrete_sofa_area(vR_init, wL_init, t_new, t_new)
gap_init = a_init - ROMIK
print(f"\n--- N={N} ({4*N} variables) ---", flush=True)
print(f"Discrete area:     {a_init:.8f} (Gap: {gap_init:+.8f})", flush=True)
red = (2.21953167 - ROMIK) / gap_init
print(f"Gap reduction vs 2.2195: {red:.1f}x", flush=True)
print(f"Elapsed:                 {time.time()-t0:.1f}s", flush=True)
