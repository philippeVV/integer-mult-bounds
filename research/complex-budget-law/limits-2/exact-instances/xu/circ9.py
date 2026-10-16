"""expXU: the exact helper-free q = 9 Hermitian (Braverman-He) schedule of expXQ as a lane-free contract circuit.

Loads ../xq/schedule9.json (order, 910 oriented arcs, R) and the F9 local codes of ../xq/code.py (p = 3).
Sessions k = 0..165 <-> back arcs R[k] = (u -> w): source port at the head w, sink port at the tail u.
F9 exactness (Phi = I) means rho = id in session indexing, with NO C/J lanes (W = 166 = S).

Provides
  * F9 local gate maps (in-arcs -> out-arcs), exact F9 transfer (must be the identity),
  * roles (a transition system: pair in-arc i with out-arc i at every gate), contract gate list, role rho,
  * T_G over GF(2^16) for random char-2 gates, and T_G mod p for random char-0 gates, on the same incidence.
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
XQ = os.path.join(HERE, '..', 'xq')   # expXQ copy (schedule9.json, code.py)
_argv = sys.argv; sys.argv = [_argv[0], '3']; sys.path.insert(0, XQ)
import code as C                                   # expXQ/code.py (verbatim expXJ field + Hermitian codes), p = 3
sys.path.pop(0); sys.argv = _argv

J = json.load(open(os.path.join(XQ, 'schedule9.json')))
q, N, k = C.q, C.N, C.k
assert q == 9 and J['q'] == 9 and list(J['D']) == list(C.D)
ORDER = J['order']; ARCS = [tuple(a) for a in J['arcs']]; R = list(J['R']); V = 2 * N; E = len(ARCS)
RSET = set(R); RIDX = {e: t for t, e in enumerate(R)}; W = len(R)
POS = {v: t for t, v in enumerate(ORDER)}
SRC = [ARCS[e][1] for e in R]                      # session enters at the head of its back arc
SNK = [ARCS[e][0] for e in R]                      # session leaves at the tail

# ---- incidence in the expXQ edge numbering: e = i*(q+1)+di joins point i and line (i+D[di]) % N (vertex N+l)
INC = [[] for _ in range(V)]; ENDS = []
for i in range(N):
    for di, dd in enumerate(C.D):
        l = (i + dd) % N; e = len(ENDS); ENDS.append((i, N + l)); INC[i].append(e); INC[N + l].append(e)
assert all(tuple(sorted(a)) == ENDS[e] for e, a in enumerate(ARCS)), 'arc e must be edge e'
INS = [[e for e in INC[v] if ARCS[e][1] == v] for v in range(V)]
OUTS = [[e for e in INC[v] if ARCS[e][0] == v] for v in range(V)]
assert all(len(INS[v]) == len(OUTS[v]) == k for v in range(V))
assert all(POS[ARCS[e][0]] < POS[ARCS[e][1]] for e in range(E) if e not in RSET), 'Gamma - R not forward in sigma'
assert all(POS[ARCS[e][0]] > POS[ARCS[e][1]] for e in R)

# ---- F9 as indices 0..8 (as in verify.py)
FQ = sorted(C.Fq); IX = {a: i for i, a in enumerate(FQ)}
ADD = np.array([[IX[int(C.add(a, b))] for b in FQ] for a in FQ], dtype=np.int64)
MUL = np.array([[IX[int(C.mul(a, b))] for b in FQ] for a in FQ], dtype=np.int64)
NEG = np.array([IX[int(C.neg(a))] for a in FQ], dtype=np.int64)
INV9 = np.zeros(q, dtype=np.int64)
for i, a in enumerate(FQ):
    if a: INV9[i] = IX[C.inv(a)]
Z9, O9 = IX[0], IX[1]
th = IX[C.gp(C.theta_star_exp)]
Gi = np.vectorize(lambda a: IX[int(a)])(C.G); Hi = np.vectorize(lambda a: IX[int(a)])(C.H)
ROW = {}                                            # ROW[(v,e)]: F9^5 generator row of edge e at vertex v
for i in range(N):
    for di, dd in enumerate(C.D):
        l = (i + dd) % N; e = i * (q + 1) + di
        ROW[(i, e)] = Gi[di]; ROW[(N + l, e)] = MUL[th if i + dd >= N else O9, Hi[di]]

# F3-coordinates of F9 (a genuine F3-basis {1, w}) -> used by label families
_F3 = (0, 1, int(C.neg(1)))
_w = next(x for x in FQ if x not in _F3)
F3COORD = {}
for a, xa in enumerate(_F3):
    for b, xb in enumerate(_F3):
        F3COORD[IX[int(C.add(xa, int(C.mul(xb, _w))))]] = (a, b)
assert len(F3COORD) == 9


def mat_inv9(A):
    n = A.shape[0]; M = np.concatenate([A.copy(), np.full((n, n), Z9, dtype=np.int64)], 1)
    for i in range(n): M[i, n + i] = O9
    for c in range(n):
        nz = [r for r in range(c, n) if M[r, c] != Z9]
        if not nz: return None
        r = nz[0]; M[[c, r]] = M[[r, c]]; M[c] = MUL[INV9[M[c, c]], M[c]]
        for r in range(n):
            if r != c and M[r, c] != Z9: M[r] = ADD[M[r], MUL[NEG[M[r, c]], M[c]]]
    return M[:, n:]


def mat_mul9(A, B):
    out = np.full((A.shape[0], B.shape[1]), Z9, dtype=np.int64)
    for t in range(A.shape[1]): out = ADD[out, MUL[A[:, t][:, None], B[t][None, :]]]
    return out


def gate9(v):
    """F9 local map: out-arc values = T_v @ in-arc values (OUTS[v] x INS[v])."""
    MI = np.array([ROW[(v, e)] for e in INS[v]]); MO = np.array([ROW[(v, e)] for e in OUTS[v]])
    Ai = mat_inv9(MI); assert Ai is not None, f'gate {v} not bijective'
    return mat_mul9(MO, Ai)


GATE9 = {v: gate9(v) for v in range(V)}


def run(gates, zero, add, mul_mat, unit):
    """Generic transfer: T[j,k] = value at sink of session j for unit input at source of session k."""
    val = {}; T = None
    rows = {}
    for v in ORDER:
        X = []
        for e in INS[v]:
            if e in RSET:
                x = np.full(W, zero, dtype=np.int64); x[RIDX[e]] = unit; X.append(x)
            else: X.append(val.pop(e))
        Y = mul_mat(gates[v], np.array(X))
        for r, e in enumerate(OUTS[v]):
            if e in RSET: rows[RIDX[e]] = Y[r]
            else: val[e] = Y[r]
    assert not val and len(rows) == W
    return np.array([rows[j] for j in range(W)])


def transfer9():
    return run(GATE9, Z9, None, mat_mul9, O9)


# ---- roles: transition system pairing INS[v][i] -> OUTS[v][i]; a role = maximal trail from a session source to a sink
def roles():
    nxt = {}
    for v in range(V):
        for a, b in zip(INS[v], OUTS[v]): nxt[a] = b
    paths = []
    for e in R:                                       # role w starts at head of R[w] (the back arc's input half)
        first = nxt[e]; path = []; cur = e
        while True:
            cur = nxt[cur]
            if cur in RSET: end = cur; break
            path.append(cur)
        paths.append((e, path, end))
    return nxt, paths


# ---------------- GF(2^16) ----------------
POLY, Q2 = 0x1100B, 1 << 16
EXP2 = np.zeros(2 * Q2, dtype=np.int64); LOG2 = np.zeros(Q2, dtype=np.int64)
_x = 1
for _i in range(Q2 - 1):
    EXP2[_i] = _x; LOG2[_x] = _i; _x <<= 1
    if _x & Q2: _x ^= POLY
EXP2[Q2 - 1:2 * Q2 - 2] = EXP2[:Q2 - 1]


def gmul(a, b):
    a = np.asarray(a, dtype=np.int64); b = np.asarray(b, dtype=np.int64)
    r = EXP2[(LOG2[a] + LOG2[b]) % (Q2 - 1)]
    return np.where((a == 0) | (b == 0), 0, r)


def ginv(a): return int(EXP2[(Q2 - 1 - LOG2[a]) % (Q2 - 1)])


def gmatmul(A, B):
    out = np.zeros((A.shape[0], B.shape[1]), dtype=np.int64)
    for t in range(A.shape[1]): out ^= gmul(A[:, t][:, None], B[t][None, :])
    return out


def grank(M):
    M = M.copy(); n, c = M.shape; r = 0
    for col in range(c):
        if r == n: break
        nz = np.nonzero(M[r:, col])[0]
        if len(nz) == 0: continue
        piv = r + nz[0]
        if piv != r: M[[r, piv]] = M[[piv, r]]
        M[r] = gmul(M[r], ginv(int(M[r, col])))
        rows = np.nonzero(M[:, col])[0]; rows = rows[rows != r]
        if len(rows): M[rows] ^= gmul(M[rows, col][:, None], M[r][None, :])
        r += 1
    return r


def transfer2(rng):
    gates = {v: rng.integers(0, Q2, size=(k, k), dtype=np.int64) for v in range(V)}
    return run(gates, 0, None, gmatmul, 1)


def Zmat2(T, betas):
    """Z[j,k] = T[j,k] (I + beta_j + beta_k) over GF(2^16) (rho = id in session indexing)."""
    m = betas[0].shape[0]; I = np.eye(m, dtype=np.int64); Zm = np.zeros((W * m, W * m), dtype=np.int64)
    for j in range(W):
        for kk in np.nonzero(T[j])[0]:
            Zm[j*m:(j+1)*m, kk*m:(kk+1)*m] = ((I + betas[j] + betas[kk]) % 2) * T[j, kk]
    return Zm


# ---------------- char 0 (mod p) ----------------
P = 33554393


def transfer_p(rng, p=P):
    def mm(G, X):
        out = np.zeros((G.shape[0], X.shape[1]), dtype=np.int64)
        for t in range(G.shape[1]): out = (out + np.outer(G[:, t], X[t])) % p
        return out
    gates = {v: rng.integers(0, p, size=(k, k), dtype=np.int64) for v in range(V)}
    return run(gates, 0, None, mm, 1)


if __name__ == '__main__':
    T9 = transfer9()
    ident = all(T9[j, kk] == (O9 if j == kk else Z9) for j in range(W) for kk in range(W))
    print(f"W = {W} sessions, V = {V} gates [10,5] over F9, E = {E} arcs, |R| = {len(R)}; "
          f"exact F9 transfer sink x source == I (rho = id on sessions, no lanes): {ident}")
