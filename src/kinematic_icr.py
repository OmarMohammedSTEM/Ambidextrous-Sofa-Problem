"""
Kinematic Instantaneous Center of Rotation (ICR) Engine for the Moving Sofa.

Aronhold-Kennedy Kinematic Theorem:
For any planar rigid body sliding against multiple fixed walls:
At every instant t, the Instantaneous Center of Rotation (ICR) c(t) must lie
on the normal line of EVERY active contact point:
    (p_i - c(t)) . t_i = 0  ==>  c(t) in line(p_i, n_i).

In a unit L-corridor:
Wall 1 (x = 0, y <= 0): Normal is (1, 0)   ==> ICR has y-coordinate y_1.
Wall 2 (y = 0, x <= 0): Normal is (0, 1)   ==> ICR has x-coordinate x_2.
Wall 3 (x = 1, y <= 1): Normal is (-1, 0)  ==> ICR has y-coordinate y_3.
Wall 4 (y = 1, x <= 1): Normal is (0, -1)  ==> ICR has x-coordinate x_4.
Inner Corner (0, 0): Normal is in cone      ==> ICR line passes through (0, 0).

Therefore:
Whenever the sofa is in 3-point contact, the ICR c(t) is COMPLETELY DETERMINED
by the contact positions without any variational search!
"""

import numpy as np

print("=== Kinematic ICR Constraint Analysis ===")

# Let's test the 5 contact states:
# State 1: Inner corner (0,0) + Outer Wall 1 (x=1) + Outer Wall 2 (y=1)
# Contact at Outer Wall 1 at y = y_o1 ==> ICR y = y_o1.
# Contact at Outer Wall 2 at x = x_o2 ==> ICR x = x_o2.
# So ICR = (x_o2, y_o1)!
# Furthermore, contact at inner corner (0,0) with sofa rotated by t:
# If corner contacts at angle t, line from (0,0) to ICR must align with the contact normal.
# This means:
# tan(t) = y_o1 / x_o2  ==>  y_o1 = x_o2 * tan(t)!

print("State 1 (Corner + Both Outer Walls):")
print("  ICR coordinates: x(t) = x_o2(t),  y(t) = y_o1(t)")
print("  Coupling relation: y(t) = x(t) * tan(t)")
print("  ==> This reduces the 2D path (xi, eta) to a SINGLE degree of freedom x(t)!")

# State 2: Inner corner (0,0) + Inner Wall 1 (x=0) + Outer Wall 1 (x=1)
# Contact at Inner Wall 1 at y = y_i1 ==> ICR y = y_i1.
# Contact at Outer Wall 1 at y = y_o1 ==> ICR y = y_o1.
# But sofa width is 1! Since x=0 and x=1 are parallel at distance 1:
# The sofa is wedged between x=0 and x=1, so y_i1 = y_o1 = y_c(t)!
# And contact with corner (0,0) gives:
# ICR x = 0 (since it touches corner at x=0)!
print("\nState 2 (Parallel corridor arm):")
print("  ICR is pinned directly to the wall axis: x(t) = 0, y(t) = y_c(t)!")

# Let's check Romik's path vs this exact kinematic prediction:
import json
d = json.load(open('C:/New folder (31)/anc/ambi_results.json'))
coef = np.array(d['B_K5_newton']['coef'])
t_eval = np.linspace(0, np.pi/2, 20)

import sys
sys.path.insert(0, 'C:/New folder (31)/anc')
import ambi
xi, eta = ambi.make_path(coef, t_eval)

print("\n--- Empirical ICR Trajectory from B_K5_newton ---")
for i in [0, 5, 9, 14, 19]:
    t = t_eval[i]
    print(f"t = {t*180/np.pi:5.1f} deg ({t:.3f} rad):  ICR = ({xi[i]:.4f}, {eta[i]:.4f})  ratio y/x = {eta[i]/max(abs(xi[i]), 1e-6):.4f}")
