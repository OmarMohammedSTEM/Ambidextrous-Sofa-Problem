import sys, time, json, numpy as np
from scipy.optimize import minimize
sys.path.insert(0, 'C:/New folder (31)/anc')
import ambi

d = json.load(open('C:/New folder (31)/anc/ambi_results.json'))
coef_sym = np.array(d['B_K5_newton']['coef'][:2+4*3]) # K=3
nz = len(coef_sym)

print(f"=== Testing Symmetry Emergence from Random Asymmetric Initialization ===")
print(f"Base symmetric area (K=3, Richardson): {ambi.fine_area(coef_sym, coef_sym)[0]:.6f}")

rng = np.random.default_rng(42)

for trial in range(3):
    noise_R = rng.normal(0, 0.05, nz)
    noise_L = rng.normal(0, 0.05, nz)
    z_init = np.concatenate([coef_sym + noise_R, coef_sym + noise_L])
    
    init_asym_norm = np.linalg.norm(z_init[:nz] - z_init[nz:])
    a_init = ambi.area(z_init[:nz], z_init[nz:], 61)
    print(f"\nTrial {trial+1}:")
    print(f"  Initial asymmetry ||v - w||: {init_asym_norm:.5f}")
    print(f"  Initial raw area:           {a_init:.6f}")
    
    # Run optimization from this asymmetric point
    f = lambda z: -ambi.area(z[:nz], z[nz:], 61)
    t0 = time.time()
    res = minimize(f, z_init, method='Powell', options=dict(maxfev=1500, xtol=1e-5, ftol=1e-8))
    
    final_R = res.x[:nz]
    final_L = res.x[nz:]
    final_asym_norm = np.linalg.norm(final_R - final_L)
    ext, a1, a2 = ambi.fine_area(final_R, final_L)
    
    print(f"  Optimized ||v - w||:         {final_asym_norm:.5f} (collapsed by {init_asym_norm / max(final_asym_norm, 1e-12):.1f}x)")
    print(f"  Optimized Richardson area:   {ext:.6f} (Romik: {ambi.ROMIK:.6f})")
    print(f"  Elapsed:                     {time.time()-t0:.1f}s")
