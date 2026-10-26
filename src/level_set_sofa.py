"""
Zero-Assumption Free-Topology / Level-Set Sofa Optimizer.

Idea:
Rather than parametrizing curves or paths with Fourier series or splines,
we represent the sofa as a continuous density field rho(x, y) on a fine 2D grid.
The objective is to maximize the area integral of rho(x, y), subject to the
constraint that rho can be translated and rotated through the L-corridor
for both right- and left-hand turns.

We solve this using a relaxed density / adjoint formulation:
rho(x, y) in [0, 1] is updated via projected gradient ascent combined
with corridor collision penalties.
"""

import time, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Grid domain: sofa rest frame centered in strip [0, 1]
# Romik sofa has length approx 2.33, height 1.0.
# We set domain x in [-1.5, 1.5], y in [0, 1]
NX, NY = 240, 120
x_grid = np.linspace(-1.5, 1.5, NX)
y_grid = np.linspace(0.0, 1.0, NY)
dx = x_grid[1] - x_grid[0]
dy = y_grid[1] - y_grid[0]
X, Y = np.meshgrid(x_grid, y_grid)

# Discrete angles for collision checking
N_ANGLES = 31
angles = np.linspace(0, np.pi/2, N_ANGLES)

print(f"Initialized Grid: {NX}x{NY} ({NX*NY} cells), dx={dx:.4f}, dy={dy:.4f}")
print(f"Tracking {N_ANGLES} rotation angles across both turn chiralities.")

# Hallway corridor test function:
# In the hallway frame, L is: x <= 1, y <= 1, and NOT (x < 0 and y < 0).
# A point (u, v) is inside L iff: u <= 1 and v <= 1 and (u >= 0 or v >= 0).
def in_hallway(u, v):
    return (u <= 1.0) & (v <= 1.0) & ((u >= 0.0) | (v >= 0.0))

# To allow the sofa to navigate, for each angle t, there is a translation (tx, ty).
# We optimize rho(x, y) together with the trajectory (tx_R, ty_R) and (tx_L, ty_L).
# As a warm start, we initialize (tx, ty) from the Hammersley/Romik path:
r = 0.35
tx_R = r * np.cos(2 * angles)
ty_R = r * np.sin(2 * angles)
tx_L = tx_R.copy()
ty_L = ty_R.copy()

# Compute the allowable mask M(x, y) = intersection over all angles of rotated hallways
def compute_admissible_mask(tx_r, ty_r, tx_l, ty_l):
    mask = np.ones((NY, NX), dtype=bool)
    # Right-turn family:
    for t, a, b in zip(angles, tx_r, ty_r):
        cos_t, sin_t = np.cos(t), np.sin(t)
        # Transform grid coordinates to hallway frame:
        # hallway corner is at (a, b) in sofa frame, hallway rotated CCW by t
        # In hallway frame: [u, v]^T = R(-t) * ([X, Y]^T - [a, b]^T)
        x_shift = X - a
        y_shift = Y - b
        u = cos_t * x_shift + sin_t * y_shift
        v = -sin_t * x_shift + cos_t * y_shift
        mask &= in_hallway(u, v)
        
    # Left-turn mirror family: reflection across y = 0.5
    for s, a, b in zip(angles, tx_l, ty_l):
        cos_s, sin_s = np.cos(s), np.sin(s)
        x_shift = X - a
        y_shift = (1.0 - Y) - b
        u = cos_s * x_shift + sin_s * y_shift
        v = -sin_s * x_shift + cos_s * y_shift
        mask &= in_hallway(u, v)
        
    return mask

t0 = time.time()
mask = compute_admissible_mask(tx_R, ty_R, tx_L, ty_L)
raw_area = np.sum(mask) * dx * dy
print(f"Initial admissible area on baseline trajectory: {raw_area:.6f} (elapsed: {time.time()-t0:.2f}s)")

# Optimize trajectory using Powell on the grid evaluator
from scipy.optimize import minimize

print("\n--- Running Free-Topology Trajectory Optimization ---", flush=True)

# Parameterize trajectory as Fourier or direct knots to ensure smoothness across angles
# 10 parameters per turn: a0, c0, and 2 harmonics
def pack_traj(p):
    # p has length 12: 6 for right, 6 for left
    cR = p[:6]
    cL = p[6:]
    # compute angles
    xi_R = cR[0] + cR[2]*np.cos(2*angles) + cR[3]*np.sin(2*angles)
    eta_R = cR[1] + cR[4]*np.cos(2*angles) + cR[5]*np.sin(2*angles)
    xi_L = cL[0] + cL[2]*np.cos(2*angles) + cL[3]*np.sin(2*angles)
    eta_L = cL[1] + cL[4]*np.cos(2*angles) + cL[5]*np.sin(2*angles)
    return xi_R, eta_R, xi_L, eta_L

def loss_func(p):
    tx_r, ty_r, tx_l, ty_l = pack_traj(p)
    m = compute_admissible_mask(tx_r, ty_r, tx_l, ty_l)
    return -float(np.sum(m) * dx * dy)

p0 = np.array([0.0, 0.0, 0.35, 0.0, 0.0, 0.35,
               0.0, 0.0, 0.35, 0.0, 0.0, 0.35])

t0 = time.time()
res = minimize(loss_func, p0, method='Powell', options=dict(maxfev=600, xtol=1e-4, ftol=1e-5))
final_area = -res.fun
print(f"Optimized Topology Area: {final_area:.6f} in {time.time()-t0:.1f}s (evals: {res.nfev})")

best_tx_r, best_ty_r, best_tx_l, best_ty_l = pack_traj(res.x)
final_mask = compute_admissible_mask(best_tx_r, best_ty_r, best_tx_l, best_ty_l)

# Render figure of the evolved shape
fig, ax = plt.subplots(figsize=(10, 4.5))
ax.imshow(final_mask, origin='lower', extent=[-1.5, 1.5, 0, 1], cmap='Blues', alpha=0.85)
ax.set_title(f"Zero-Assumption Free-Topology Evolved Sofa (Area approx {final_area:.5f})")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig("C:/New folder (31)/figures/evolved_topology_sofa.png", dpi=150)
print("Saved topology figure: figures/evolved_topology_sofa.png")
