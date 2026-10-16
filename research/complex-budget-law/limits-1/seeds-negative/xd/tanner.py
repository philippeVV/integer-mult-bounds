"""expXD: Tanner codes on H(q) (or H(q)-M) with square [t,t/2] local codes.
Edge values x_{p,l}. Point p: x_{p,.} = c_p^T G_p (k x deg). Line l: H_l x_{.,l} = 0.
dim W = nullity of the line-check matrix in the NP*k unknowns c_p.
Natural (geometric) gauge: point side coords = coordinates of the line's normalized Plucker vector in a basis of
p ^ pi_p; line side coords = coordinates of the point's normalized vector in a basis of l.  Both are G2(q)-equivariant up to
GL2 at each vertex, so dim W is gauge-invariant w.r.t. the choice of normalized reps."""
import sys, itertools, random, collections, math, numpy as np, numba as nb
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from hexgeo import build, incidence, graph, girth, GF

@nb.njit(cache=True)
def _rref(A, add, mul, neg, inv, piv):
    nr, nc = A.shape; r = 0
    for c in range(nc):
        if r >= nr: break
        p = -1
        for i in range(r, nr):
            if A[i, c] != 0: p = i; break
        if p < 0: continue
        if p != r:
            for k in range(nc):
                t = A[r, k]; A[r, k] = A[p, k]; A[p, k] = t
        iv = inv[A[r, c]]
        for k in range(c, nc): A[r, k] = mul[iv, A[r, k]]
        for i in range(nr):
            if i != r:
                f = A[i, c]
                if f != 0:
                    nf = neg[f]
                    for k in range(c, nc):
                        if A[r, k] != 0: A[i, k] = add[A[i, k], mul[nf, A[r, k]]]
        piv[r] = c; r += 1
    return r

def rref(F, M):
    A = np.ascontiguousarray(np.array(M, dtype=np.int16)); piv = np.zeros(A.shape[1], dtype=np.int64)
    r = _rref(A, F.add, F.mul, F.neg, F.inv, piv); return A[:r], list(piv[:r])
def rank(F, M): return len(rref(F, M)[1])
def nullspace(F, M):
    R, piv = rref(F, M); nc = M.shape[1]; ps = set(piv); free = [c for c in range(nc) if c not in ps]
    B = np.zeros((len(free), nc), dtype=np.int16)
    for t, f in enumerate(free):
        B[t, f] = 1
        for i, pc in enumerate(piv): B[t, pc] = F.neg[R[i, f]]
    return B

def coords(Fq, vecs, target):
    """find (a,b) in Fq^2 with a*vecs[0]+b*vecs[1] == target (vectors as tuples)"""
    A, M = Fq.add, Fq.mul; q = Fq.q
    for a in range(q):
        for b in range(q):
            if all(int(A[M[a, x], M[b, y]]) == z for x, y, z in zip(vecs[0], vecs[1], target)): return a, b
    raise ValueError

def geo_coords(H):
    """returns ptc[(p,l)] = coords of line l at point p ; lnc[(p,l)] = coords of point p on line l (over F_q)"""
    Fq = H['F']; NP, NL = incidence(H)
    lnc = {}
    for l, mem in enumerate(H['lines']):
        u, w = H['pts'][mem[0]], H['pts'][mem[1]]
        for p in mem: lnc[(p, l)] = coords(Fq, (u, w), H['pts'][p])
    ptc = {}
    for p in range(NP):
        ls = H['adjP'][p]; u, w = H['plk'][ls[0]], H['plk'][ls[1]]
        for l in ls: ptc[(p, l)] = coords(Fq, (u, w), H['plk'][l])
    return ptc, lnc

def embed(Fq, K):
    """field embedding Fq -> K as int arrays (Fq, K both from gf.GF with same p)"""
    if Fq.q == Fq.p: return np.arange(Fq.q)
    # find generator image: Fq built with irreducible; find root in K of its minimal poly by brute force over homs
    for g in range(K.q):
        img = np.zeros(Fq.q, dtype=np.int64)
        # basis 1, x, ..., x^{e-1}: element a = sum d_i x^i
        pows = [1]
        for i in range(1, Fq.e): pows.append(int(K.mul[pows[-1], g]))
        for a in range(Fq.q):
            v = 0; t = a
            for i in range(Fq.e):
                d = t % Fq.p; t //= Fq.p
                for _ in range(d): v = int(K.add[v, pows[i]])
            img[a] = v
        if len(set(img)) == Fq.q and all(img[int(Fq.mul[a, b])] == K.mul[img[a], img[b]] for a in range(Fq.q) for b in range(Fq.q)):
            return img
    raise ValueError

def randmds(K, k, n, rng):
    while True:
        G = np.array([[rng.randrange(K.q) for _ in range(n)] for _ in range(k)], dtype=np.int16)
        if all(rank(K, G[:, list(c)]) == k for c in itertools.combinations(range(n), k)): return G

def frob(K, x, j):
    for _ in range(j): x = K.pw(x, K.p)
    return x

class Tanner:
    def __init__(s, H, K, drop=(), sp=0, sl=0, mode='nat', seed=0, dualP=False, dualL=False, kP=None, kL=None):
        s.H, s.K = H, K; Fq = H['F']; NP, NL = incidence(H); s.NP, s.NL = NP, NL
        emb = embed(Fq, K); rng = random.Random(seed); s.drop = set(drop)
        ptc, lnc = geo_coords(H)
        s.pcols = [[l for l in H['adjP'][p] if (p, l) not in s.drop] for p in range(NP)]
        s.lcols = [[p for p in H['lines'][l] if (p, l) not in s.drop] for l in range(NL)]
        deg = len(s.pcols[0]); assert all(len(c) == deg for c in s.pcols + s.lcols); s.deg = deg
        kP = kP or deg // 2; kL = kL or deg // 2; s.kP, s.kL = kP, kL; s.k = kP
        assert kP + kL == deg
        def homcol(ab, j, kk):
            a, b = int(emb[ab[0]]), int(emb[ab[1]]); a, b = frob(K, a, j), frob(K, b, j)
            return [int(K.mul[K.pw(a, kk - 1 - i), K.pw(b, i)]) for i in range(kk)]
        def rnz(): return rng.randrange(1, K.q)
        s.GP = []
        for p in range(NP):
            if mode == 'randgate':
                G = randmds(K, kP, deg, rng)
            else:
                G = np.array([homcol(ptc[(p, l)], sp, kP) for l in s.pcols[p]], dtype=np.int16).T.copy()
                if mode == 'randgauge':
                    for c in range(deg):
                        g = rnz(); G[:, c] = K.mul[g, G[:, c]]
            if dualP: G = nullspace(K, G)
            s.GP.append(G)
        s.HL = []
        for l in range(NL):
            if mode == 'randgate':
                G = randmds(K, kL, deg, rng)
            else:
                G = np.array([homcol(lnc[(p, l)], sl, kL) for p in s.lcols[l]], dtype=np.int16).T.copy()
            if dualL: G = nullspace(K, G)
            s.HL.append(nullspace(K, G))
        s.E = NP * deg
        s.pidx = {(p, l): c for p in range(NP) for c, l in enumerate(s.pcols[p])}
    def mds_ok(s):
        K = s.K
        def mds(G):
            kk, n = G.shape
            return all(rank(K, G[:, list(c)]) == kk for c in itertools.combinations(range(n), kk))
        return all(mds(G) for G in s.GP[:20]) and all(mds(nullspace(K, Hh)) for Hh in s.HL[:20])
    def check_matrix(s):
        K = s.K; k = s.kP
        rows = np.zeros((s.NL * (s.deg - s.kL), s.NP * k), dtype=np.int16); r = 0
        for l in range(s.NL):
            for h in s.HL[l]:
                for t, p in enumerate(s.lcols[l]):
                    if h[t] == 0: continue
                    col = s.GP[p][:, s.pidx[(p, l)]]
                    rows[r, p * k:(p + 1) * k] = K.add[rows[r, p * k:(p + 1) * k], K.mul[h[t], col]]
                r += 1
        return rows
    def dimW(s):
        s.Kern = nullspace(s.K, s.check_matrix()); return s.Kern.shape[0]

def a_of(S, E, b=11):
    if S <= 0: return float('-inf')
    return 1 - math.log((E - S) / S) / math.log(b)

def matching(H, rng):
    """random perfect matching in point-line incidence (regular bipartite) via randomized augmenting paths"""
    NP, NL = incidence(H)
    adj = [list(a) for a in H['adjP']]
    for a in adj: rng.shuffle(a)
    mL = [-1] * NL; mP = [-1] * NP
    order = list(range(NP)); rng.shuffle(order)
    sys.setrecursionlimit(100000)
    for p0 in order:
        # BFS augment
        prev = {}; dq = collections.deque([p0]); seen = {p0}; found = -1
        while dq and found < 0:
            p = dq.popleft()
            for l in adj[p]:
                if l in prev: continue
                prev[l] = p
                if mL[l] < 0: found = l; break
                if mL[l] not in seen: seen.add(mL[l]); dq.append(mL[l])
        assert found >= 0
        l = found
        while True:
            p = prev[l]; nxt = mP[p]; mP[p] = l; mL[l] = p
            if p == p0: break
            l = nxt
    return [(p, mP[p]) for p in range(NP)]
