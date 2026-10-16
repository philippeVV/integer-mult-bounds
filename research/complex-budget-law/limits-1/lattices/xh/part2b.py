"""expXH part 2b: Lorentzian Hermitian 'E_n' port systems over O = Z[i] or Z[omega] (u = 2).
Ports: v in O^{n,1}, H(v,v) = -|v0|^2 + sum |vj|^2 = 2, optionally H(v,k) = 0 with k = (3;1^n).  Heights |v0|^2 <= hmax.
Valid: every nonzero H_ab in 1 + 2O (labels (v,1), gauge absorbed by allowing all unit multiples as candidates).
Reports: #candidates, #candidates that are unit multiples of real vectors, max valid set V (CP-SAT), complex rank,
h_R, h_C (= 2 h_R by the tensor law), gamma_C of the complex coordinate decoder (n+1 terms, cost = 2*rank B_R[S,S]),
and the real I_{n,1} control (same constraints, O = Z).
usage: part2b.py ring(gauss|eis|real) n hmax [k]"""
import sys, itertools
import numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from lib import rank_f2_rows, bits_of, max_clique

ring, n, hmax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
usek = len(sys.argv) > 4 and sys.argv[4] == 'k'
om = complex(-0.5, 3 ** 0.5 / 2)
if ring == 'gauss':
    units = [1, -1, 1j, -1j]
    elems = [complex(a, b) for a in range(-2, 3) for b in range(-2, 3)]
elif ring == 'eis':
    units = [1, -1, om, -om, om * om, -om * om]
    elems = [a + b * om for a in range(-2, 3) for b in range(-2, 3)]
else:
    units = [1, -1]
    elems = [complex(a, 0) for a in range(-2, 3)]
nrm = lambda z: round(abs(z) ** 2, 9)
elems = [z for z in elems if nrm(z) <= 4]
elems = list({(round(z.real, 9), round(z.imag, 9)): z for z in elems}.values())
cands = []
for v0 in elems:
    t = 2 + nrm(v0)
    if nrm(v0) > hmax: continue
    cur = []
    def rec(i, q, s):
        if i == n:
            if abs(q - t) < 1e-9 and (not usek or abs(s - 3 * v0) < 1e-9):
                cands.append([v0] + cur.copy())
            return
        for z in elems:
            if q + nrm(z) <= t + 1e-9:
                cur.append(z); rec(i + 1, q + nrm(z), s + z); cur.pop()
    rec(0, 0.0, 0)
# one representative per {+-1} (sign does not change label); keep other unit multiples as separate candidates
reps = {}
for v in cands:
    key = tuple((round(z.real, 6), round(z.imag, 6)) for z in v)
    neg = tuple((round(-z.real, 6), round(-z.imag, 6)) for z in v)
    if neg not in reps: reps[key] = v
X = np.array(list(reps.values()), dtype=complex)
eps = np.array([-1] + [1] * n)
def is_real_mult(v):
    for u in units:
        w = v / u
        if np.allclose(w.imag, 0): return True
    return False
nreal = sum(is_real_mult(v) for v in X)
H = (X * eps) @ X.conj().T
assert np.allclose(np.diag(H), 2)
def in_1p2O(z):
    if ring == 'gauss':
        a, b = round(z.real), round(z.imag); return a % 2 == 1 and b % 2 == 0
    if ring == 'eis':
        q = 2 * z.imag / 3 ** 0.5; p = z.real + q / 2; P, Q = round(p), round(q)
        assert abs(p - P) < 1e-7 and abs(q - Q) < 1e-7
        return (P - 1) % 2 == 0 and Q % 2 == 0
    return round(z.real) % 2 == 1
Z0 = np.isclose(np.abs(H), 0)
adj = Z0.copy()
V0 = len(X)
for a in range(V0):
    for b in range(a + 1, V0):
        if not Z0[a, b] and in_1p2O(H[a, b]):
            adj[a, b] = adj[b, a] = True
np.fill_diagonal(adj, False)
cl, opt, ub = max_clique(adj, seconds=120)
S = np.array(cl); V = len(S)
Hs = H[np.ix_(S, S)]
rC = np.linalg.matrix_rank(Hs, tol=1e-8)
BR = (np.isclose(np.abs(Hs), 0) | np.eye(V, dtype=bool)).astype(np.uint8)
f2 = lambda M: rank_f2_rows(bits_of(M))
hR = f2(BR)
CR = 0
for c in range(n + 1):
    Sc = np.nonzero(~np.isclose(np.abs(X[S, c]), 0))[0]
    if len(Sc): CR += f2(BR[np.ix_(Sc, Sc)])
u = 1 if ring == 'real' else 2
h = u * hR
gam = CR / V
allreal = all(is_real_mult(v) for v in X[S])
print(f"{ring:5s} n={n} hmax={hmax} k={'(3;1^n)' if usek else 'none'}: cands={V0} (unit-multiples of real: {nreal}) | "
      f"maxvalid V={V} opt={opt} ub={ub} allRealMult={allreal} rank_C={rC} h_R={hR} u={u} h={h} C={u * CR} gamma={gam:.3f} "
      f"merit gamma*(h-1)/u={gam * (h - 1) / u:.2f}", flush=True)

# limits-1 replay: with k = (3;1^n), every valid port is a unit multiple of a real E_n root and the max valid set equals
# the real count (4, 10, 20, 36, 63 for n = 3..7)
if usek and ring != 'real':
    assert allreal and opt and V == {3: 4, 4: 10, 5: 20, 6: 36, 7: 63}[n], (allreal, opt, V)
print('PASS_EXPXH_PART2B')
