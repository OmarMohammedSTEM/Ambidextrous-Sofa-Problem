import numpy as np

# Dan Romik's exact constants from romik.js & sofa-math.js
beta = 0.289653820817320941743521611736

k11 = 0.124712637587267758739932415305
k12 = 0.5
a1  = 0.875287362412732241260067584695
a2  = 0.0

k61 = -0.167049816550309655013423446260
k62 = 0.5
f1  = 1.20293890815691138907022280034
f2  = -0.498273610464875672029397859080

k51 = -0.458812270687887068766779307825
k52 = 0.5
e1  = 0.875287362412732241260067584695
e2  = 0.0

def rot(t, vx, vy):
    c, s = np.cos(t), np.sin(t)
    return c * vx - s * vy, s * vx + c * vy

def romik_path_single(t):
    if t < beta:
        vx = a1 * np.cos(t) + a2 * np.sin(t) - 1.0
        vy = -a2 * np.cos(t) + a1 * np.sin(t) - 0.5
        rx, ry = rot(t, vx, vy)
        return rx + k11, ry + k12
    elif t < np.pi/2 - beta:
        vx = f1 * np.cos(t/2) + f2 * np.sin(t/2) - 1.0
        vy = -f2 * np.cos(t/2) + f1 * np.sin(t/2) - 1.0
        rx, ry = rot(t, vx, vy)
        return rx + k61, ry + k62
    else:
        vx = e1 * np.cos(t) + e2 * np.sin(t) - 0.5
        vy = -e2 * np.cos(t) + e1 * np.sin(t) - 1.0
        rx, ry = rot(t, vx, vy)
        return rx + k51, ry + k52

def romik_path(ts):
    pts = [romik_path_single(t) for t in ts]
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])

if __name__ == '__main__':
    import sys
    sys.path.insert(0, 'C:/New folder (31)/anc')
    from shapely.geometry import Polygon
    from shapely import affinity

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
        # Mirrored path for left turn:
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

    ROMIK_EXACT = 1.644955218425440

    print("=== Testing Dan Romik's Exact Analytic Path ===")
    print(f"Romik Closed-Form Theoretical Area: {ROMIK_EXACT:.12f}")
    
    raw_areas = []
    Ns = [81, 161, 321, 641, 1281]
    for N in Ns:
        S = romik_sofa_poly(N)
        a = S.area
        raw_areas.append(a)
        xs, ys = np.array(S.exterior.xy)
        L = xs.max() - xs.min()
        print(f"N={N:4d}: Raw Discrete Area = {a:.8f} (Bias: {a - ROMIK_EXACT:+.8f}), Length = {L:.7f}")
    
    print("\n--- Richardson Extrapolation on Romik's Exact Path ---")
    ext_161 = 2 * raw_areas[1] - raw_areas[0]
    ext_321 = 2 * raw_areas[2] - raw_areas[1]
    ext_641 = 2 * raw_areas[3] - raw_areas[2]
    ext_1281 = 2 * raw_areas[4] - raw_areas[3]
    print(f"Richardson (161 vs  81):  {ext_161:.10f}  (Error: {ext_161 - ROMIK_EXACT:+.2e})")
    print(f"Richardson (321 vs 161):  {ext_321:.10f}  (Error: {ext_321 - ROMIK_EXACT:+.2e})")
    print(f"Richardson (641 vs 321):  {ext_641:.10f}  (Error: {ext_641 - ROMIK_EXACT:+.2e})")
    print(f"Richardson (1281 vs 641): {ext_1281:.10f}  (Error: {ext_1281 - ROMIK_EXACT:+.2e})")
