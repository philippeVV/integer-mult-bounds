"""|S| = 1 geometry: s = line 0. E = orth(s) (63), N = nonorth(s) (56).
For each class h: rank of R[O_h cap E] inside r_s^perp; candidate generic-line Def sets."""
import sys, numpy as np, itertools, collections
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from e8 import *
from maxdef import span2
R = np.array(ROOT, float) / 2          # true coords, r.r = 2
s = 0
E = [a for a in range(NL) if IP[a][s] == 0 and a != s]; N = [a for a in range(NL) if IP[a][s] != 0 and a != s]
CL = {h: [a for a in range(NL) if dotq(LAB[a], h)] for h in range(1, 511)}
def rk(L): return np.linalg.matrix_rank(R[L]) if L else 0
cnt = collections.Counter(); cand = {}
for h, O in CL.items():
    OE = [a for a in O if a in E]; r = rk(OE)
    cnt[(len(OE), r, s in O)] += 1
    if r <= 6: cand[h] = OE
print('(|O cap E|, rank, s in O) ->', sorted(cnt.items()))
print('labels of C(s) span dim', len(span2([LAB[a] for a in N + [s]])).bit_length() - 1)
# subspaces V_h = (span R[O cap E])^perp inside r_s^perp, as orthonormal bases
P = np.eye(8) - np.outer(R[s], R[s]) / 2
def perpbasis(L):
    M = np.vstack([R[L], R[s]]) if L else R[s][None]
    u, sv, vt = np.linalg.svd(M); r = (sv > 1e-9).sum(); return vt[r:]
V = {h: perpbasis(OE) for h, OE in cand.items()}
print('candidates', len(cand), 'dims of V_h', collections.Counter(len(v) for v in V.values()))
# closure of intersections: subspaces as frozenset of h-sets that contain them
def inter(A, B):
    M = np.vstack([A, B]); # intersection of row spaces: null of [A; -B]^T
    K = np.linalg.svd(np.vstack([A, -B]).T)[2]; sv = np.linalg.svd(np.vstack([A, -B]).T, compute_uv=False)
    r = (sv > 1e-9).sum(); null = K[r:]
    if len(null) == 0: return None
    W = null[:, :len(A)] @ A
    u, sv2, vt = np.linalg.svd(W); r2 = (sv2 > 1e-9).sum(); return vt[:r2] if r2 else None
def Dof(W):
    return frozenset(h for h, Vh in V.items() if np.linalg.matrix_rank(np.vstack([Vh, W]), 1e-8) == len(Vh))
seen = {}
frontier = [V[h] for h in V]
for W in frontier:
    D = Dof(W)
    if D not in seen: seen[D] = W
level = list(seen.items())
for it in range(4):
    new = {}
    for (D1, W1), (D2, W2) in itertools.combinations(level, 2):
        if D1 <= D2 or D2 <= D1: continue
        W = inter(W1, W2)
        if W is None: continue
        D = Dof(W)
        if D not in seen and D not in new: new[D] = W
    print('iter', it, 'new', len(new), flush=True)
    if not new: break
    seen.update(new); level = list(seen.items())
dims = collections.Counter((len(span2(sorted(D))).bit_length() - 1, len(W)) for D, W in seen.items())
print('(dim span D, dim W) counts', sorted(dims.items()))
for D, W in seen.items():
    d = len(span2(sorted(D))).bit_length() - 1
    assert len(span2(sorted(D))) - 1 == len(D), 'D not closed (not a subspace)'
mx = max(len(D) for D in seen)
print('max |D|', mx, 'dim', mx.bit_length() - 1 if mx else 0)
