"""Ambidextrous moving sofa: path-space optimization.

Sofa frame: sofa fixed in strip 0<=y<=1.
Right turn: hallway H0 (inner corner at origin, outer corner (1,1), arms -> -x,-y)
rotated CCW by t about corner, corner translated to v(t).  (Validated vs Hammersley.)
Left turn: mirror family, reflection across y=1/2, independent path w(s).
Sofa = STRIP  cap_t L_t(v)  cap_s rho(L_s(w)).  Maximize area over Fourier coeffs.
"""
import numpy as np, json, time, sys
from shapely.geometry import Polygon
from shapely import affinity
from scipy.optimize import minimize

ARM = 6.0
H0 = Polygon([(1,1),(-ARM,1),(-ARM,0),(0,0),(0,-ARM),(1,-ARM)])
STRIP = Polygon([(-ARM,0),(ARM,0),(ARM,1),(-ARM,1)])
PI2 = np.pi/2

def make_path(coef, ts):
    """coef = [a0, c0, A1..AK, B1..BK, C1..CK, D1..DK]
    xi = a0 + sum Ak cos2kt + Bk sin2kt ; eta = c0 + sum Ck cos2kt + Dk sin2kt"""
    coef = np.asarray(coef, float)
    K = (len(coef)-2)//4
    A = coef[2:2+K]; B = coef[2+K:2+2*K]; C = coef[2+2*K:2+3*K]; D = coef[2+3*K:2+4*K]
    xi  = np.full_like(ts, coef[0]); eta = np.full_like(ts, coef[1])
    for k in range(1, K+1):
        ck, sk = np.cos(2*k*ts), np.sin(2*k*ts)
        xi  += A[k-1]*ck + B[k-1]*sk
        eta += C[k-1]*ck + D[k-1]*sk
    return xi, eta

def sofa_poly(coefR, coefL, N):
    ts = np.linspace(0.0, PI2, N)
    S = STRIP
    for coef, mirror in ((coefR, False), (coefL, True)):
        xi, eta = make_path(coef, ts)
        for t, a, b in zip(ts, xi, eta):
            Ht = affinity.translate(
                affinity.rotate(H0, t, origin=(0,0), use_radians=True), a, b)
            if mirror:
                Ht = affinity.scale(Ht, 1.0, -1.0, origin=(0.0, 0.5))
            S = S.intersection(Ht)
            if S.is_empty:
                return None
    if S.geom_type != 'Polygon':
        polys = [g for g in getattr(S, 'geoms', [S]) if g.geom_type == 'Polygon']
        if not polys:
            return None
        S = max(polys, key=lambda p: p.area)
    return S

def area(coefR, coefL, N):
    S = sofa_poly(coefR, coefL, N)
    return 0.0 if S is None else S.area

def fine_area(coefR, coefL, N=241):
    """Richardson-extrapolated area (bias is ~c/N, verified on Hammersley)."""
    a1 = area(coefR, coefL, N)
    a2 = area(coefR, coefL, 2*N-1)
    return 2*a2 - a1, a1, a2

def pad(coef, K):
    """re-embed coef vector into a larger harmonic space K."""
    coef = np.asarray(coef, float)
    k0 = (len(coef)-2)//4
    out = np.zeros(2+4*K); out[:2] = coef[:2]
    for blk in range(4):
        out[2+blk*K : 2+blk*K+k0] = coef[2+blk*k0 : 2+(blk+1)*k0]
    return out

ROMIK = 1.644955218  # Romik 2017 ambidextrous sofa

if __name__ == '__main__':
    results = {}
    # ---------- Stage A: symmetric double-Hammersley baseline ----------
    t0 = time.time()
    best_r, best_a = None, -1
    for r in np.linspace(0.25, 0.75, 21):
        c = np.array([0,0, r, 0, 0, r])       # K=1: xi=r cos2t, eta=r sin2t
        a = area(c, c, 81)
        if a > best_a: best_a, best_r = a, r
    print(f'[A] double-Hammersley best: r={best_r:.4f} raw81={best_a:.6f}', flush=True)
    c = np.array([0,0, best_r, 0, 0, best_r])
    ext,_,_ = fine_area(c, c)
    print(f'[A] Richardson area = {ext:.6f}  (Romik = {ROMIK})  [{time.time()-t0:.0f}s]', flush=True)
    results['A_baseline'] = dict(r=best_r, area=ext)

    # ---------- Stage B: symmetric Fourier ladder ----------
    N_SEARCH = 61
    x = np.array([0,0, best_r, 0, 0, best_r], float)
    for K, mf in ((1, 800), (2, 1600), (3, 2600)):
        x = pad(x, K)
        t0 = time.time()
        f = lambda p: -area(p, p, N_SEARCH)
        res = minimize(f, x, method='Powell',
                       options=dict(maxfev=mf, xtol=1e-6, ftol=1e-9))
        x = res.x
        ext, a1, a2 = fine_area(x, x)
        print(f'[B] K={K}: raw{N_SEARCH}={-res.fun:.6f}  Richardson={ext:.6f} '
              f'(gap to Romik {ext-ROMIK:+.2e})  nfev={res.nfev} [{time.time()-t0:.0f}s]', flush=True)
        results[f'B_K{K}'] = dict(coef=list(x), area=ext)

    # ---------- Stage C: break the up-down (chirality) symmetry ----------
    t0 = time.time()
    K = 3
    z0 = np.concatenate([x, x])
    nz = len(x)
    f2 = lambda z: -area(z[:nz], z[nz:], N_SEARCH)
    res2 = minimize(f2, z0, method='Powell',
                    options=dict(maxfev=4000, xtol=1e-6, ftol=1e-9))
    ext2, _, _ = fine_area(res2.x[:nz], res2.x[nz:])
    ext_sym, _, _ = fine_area(x, x)
    print(f'[C] asymmetric: Richardson={ext2:.6f}  vs symmetric {ext_sym:.6f} '
          f'(gain {ext2-ext_sym:+.2e})  nfev={res2.nfev} [{time.time()-t0:.0f}s]', flush=True)
    results['C_asym'] = dict(coefR=list(res2.x[:nz]), coefL=list(res2.x[nz:]), area=ext2)

    # ---------- Stage D: targeted probes at the symmetric optimum ----------
    # notch-offset mode (a0 of the mirrored path) and a few random directions
    print('[D] directional probes (central FD, N=161, h=0.01):', flush=True)
    base = area(x, x, 161)
    dirs = {'offset a0(L)': ('L', 0)}
    for name, (side, idx) in dirs.items():
        for h in (0.01, 0.03):
            xp = x.copy(); xm = x.copy()
            xp[idx] += h; xm[idx] -= h
            ap = area(x, xp, 161) if side=='L' else area(xp, x, 161)
            am = area(x, xm, 161) if side=='L' else area(xm, x, 161)
            print(f'    {name} h={h}: dA={ (ap-am)/(2*h):+.5f}  '
                  f'curv={(ap+am-2*base)/h**2:+.3f}', flush=True)
    rng = np.random.default_rng(1)
    worst = 0.0
    for j in range(8):
        d = rng.standard_normal(2*nz); d /= np.linalg.norm(d)
        h = 0.02
        zp = z0 + h*np.concatenate([d[:nz], d[nz:]])
        ap = area(zp[:nz], zp[nz:], 161)
        print(f'    random dir {j}: A(x+h d)-A(x) = {ap-base:+.2e}', flush=True)
        worst = max(worst, ap-base)
    print(f'[D] best random perturbation gain: {worst:+.2e}', flush=True)

    with open('data/ambi_results.json','w') as fh:
        json.dump(results, fh, indent=1)
    print('DONE', flush=True)
