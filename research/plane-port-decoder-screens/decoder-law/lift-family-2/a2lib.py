"""Non-concurrent A2-triple lift, normal form (gauge c0 = 0, nu_1 = nu_2 = 0, nu_3 = 1 by GL9 shift + 9th-column scaling):
  X: x_t = (r_t, 0) (t in T = V minus S), x_s = xi_s (c_s, 1), c_si = mu_i r_si (i = 1, 2), c_s3 = mu_3 r_s3 + n,
  Y: y_t = (r_t/2, eta_t), eta_t = -(n.r_t)/2 on T_3 (lines orth to s3 only), 0 elsewhere on T;  y_s = (beta_s r_s, eta_s).
  n = the vector of P = span(r_s1, r_s2) orthogonal to r_s3 (scaled so r_s1.n = 1).
The Y_T part has NO free parameter. Free: mu_i, xi_s, and the projective directions of y_s (theta_s = eta_s/beta_s)."""
import sys, json, itertools
from fractions import Fraction as Fr
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from kfam import RT, dot, kernel, normline, CLASS
from e8 import *
orb = json.load(open('./orbits.json')); S = next(o['rep'] for o in orb['3'] if o['rankS'] == 2)
s1, s2, s3 = S
T = [a for a in range(NL) if a not in S]
def orth(a, b): return IP[a][b] == 0
T0 = [t for t in T if all(orth(t, s) for s in S)]
Ti = {s: [t for t in T if orth(t, s) and not all(orth(t, x) for x in S)] for s in S}
assert all(any(orth(t, s) for s in S) for t in T), 'F not empty'
# n in P, n . r_s3 = 0
a1, a2, a3 = RT[s1], RT[s2], RT[s3]
# n = x a1 + y a2 with n.a3 = 0
d1, d2 = dot(a1, a3), dot(a2, a3)
n = [d2 * u - d1 * v for u, v in zip(a1, a2)]
assert dot(n, a3) == 0 and any(n)
YT = {}
for t in T:
    eta = -dot(n, RT[t]) / 2 if t in Ti[s3] else Fr(0)
    YT[t] = [c / 2 for c in RT[t]] + [eta]
def ylines_T():
    """for each class h: exact kernel basis of Y[O_h cap T]."""
    out = {}
    for h, O in CLASS.items():
        out[h] = kernel([YT[a] for a in O if a in YT])
    return out
