import sys
sys.path.insert(0, 'C:/New folder (31)/anc')
import numpy as np
from shapely.geometry import Polygon
from shapely import affinity
from test_romik_exact_path import romik_path

ROMIK_EXACT = 1.644955218425440
ARM = 6.0
H0 = Polygon([(1,1),(-ARM,1),(-ARM,0),(0,0),(0,-ARM),(1,-ARM)])
STRIP = Polygon([(-ARM,0),(ARM,0),(ARM,1),(-ARM,1)])

def romik_sofa_poly(N):
    ts = np.linspace(0, np.pi/2, N)
    xi, eta = romik_path(ts)
    S = STRIP
    for t, px, py in zip(ts, xi, eta):
        Ht = affinity.translate(
            affinity.rotate(H0, float(t), origin=(0,0), use_radians=True),
            float(px), float(py))
        S = S.intersection(Ht)
        if S.is_empty: return None
    for t, px, py in zip(ts, xi, eta):
        Hs = affinity.translate(
            affinity.rotate(H0, float(t), origin=(0,0), use_radians=True),
            float(px), float(py))
        Hs = affinity.scale(Hs, 1.0, -1.0, origin=(0.0, 0.5))
        S = S.intersection(Hs)
        if S.is_empty: return None
    if S.geom_type != 'Polygon':
        polys = [g for g in getattr(S, 'geoms', [S]) if g.geom_type == 'Polygon']
        S = max(polys, key=lambda p: p.area)
    return S

if __name__ == '__main__':
    print("=== Convergence of Polygon Enclosure to Continuum Romik Constant ===")
    Ns = [50, 100, 200, 400, 800, 1600]
    areas = []
    for N in Ns:
        p = romik_sofa_poly(N)
        a = p.area
        areas.append(a)
        bias = a - ROMIK_EXACT
        c1 = bias * N
        print(f"N = {N:4d}: Area = {a:.10f} | Bias = {bias:+.8f} | c1 = N*Bias = {c1:.6f}", flush=True)

    print("\n=== Romberg Tableau on Romik Enclosure Sequence ===")
    # Successive Romberg extrapolation steps:
    # A_k^{(m)} = (2^m * A_{k+1}^{(m-1)} - A_k^{(m-1)}) / (2^m - 1)
    # We will do this for dyadic sequence: N = 100, 200, 400, 800, 1600
    sub_A = areas[1:] # 100, 200, 400, 800, 1600
    R = [sub_A]
    for m in range(1, len(sub_A)):
        prev = R[-1]
        factor = 2**m
        curr = [(factor * prev[i+1] - prev[i]) / (factor - 1.0) for i in range(len(prev)-1)]
        R.append(curr)
        print(f"Order m={m}:")
        for val in curr:
            diff = val - ROMIK_EXACT
            print(f"   Value = {val:.14f} | Error to Romik = {diff:+.2e}", flush=True)
