"""expQA core: masked-identity families and min-rank tools.

Mask: N_aa = 1, N_ab = 0 when ports a,b are NOT orthogonal (required), free when orthogonal.
free[a,b] = True on orthogonal pairs (a != b).
"""
import sys, itertools
import numpy as np
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # fam.py (expNM), minrank.py (expNQ), copied here


def family(name):
    """name: 'e8' | 'pg:n:q:dot|herm[:plus1]'. returns free (bool VxV), theta, h, label, points"""
    if name == 'e8':
        from minrank import e8_family
        free, th, h, lab = e8_family()
        return free, th, h, lab, None
    _, n, q, kind, *rest = name.split(':')
    n, q = int(n), int(q); u = q.bit_length() - 1
    from fam import build
    F, B, P = build(u, n, kind, bool(rest))
    V = len(P)
    # vectorised form evaluation
    Pm = np.array(P, dtype=np.int64)
    mul = np.array([[F.mul(a, b) for b in range(F.q)] for a in range(F.q)], dtype=np.int64)
    if kind == 'herm':
        Qr = 1 << (u // 2)
        conj = np.array([F.pw(b, Qr) for b in range(F.q)], dtype=np.int64)
        Pc = conj[Pm]
    else:
        Pc = Pm
    Bm = np.zeros((V, V), dtype=np.int64)
    for i in range(n):
        Bm ^= mul[Pm[:, i][:, None], Pc[:, i][None, :]]
    free = (Bm == 0)
    np.fill_diagonal(free, False)
    return free, 1.0 / n, u * n, f"PG({n-1},{q}) {kind}{' +1' if rest else ''}", P


def frames(free, n):
    """all free n-cliques (frames)."""
    V = free.shape[0]
    adj = [set(np.nonzero(free[a])[0].tolist()) for a in range(V)]
    cl = []
    def ext(c, cand):
        if len(c) == n: cl.append(tuple(c)); return
        for b in sorted(cand):
            if b > c[-1]: ext(c + [b], cand & adj[b])
    for a in range(V): ext([a], adj[a])
    return cl


def clique_partition(free, tl=60, workers=2, kmax=None):
    """min clique partition of the free graph via CP-SAT over all cliques of size <= kmax (default omega guess)."""
    from ortools.sat.python import cp_model
    V = free.shape[0]
    adj = [set(np.nonzero(free[a])[0].tolist()) for a in range(V)]
    cl = []
    def ext(c, cand):
        cl.append(tuple(c))
        if kmax and len(c) == kmax: return
        for b in sorted(cand):
            if b > c[-1]: ext(c + [b], cand & adj[b])
    for a in range(V): ext([a], adj[a])
    m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in cl]
    inc = [[] for _ in range(V)]
    for i, c in enumerate(cl):
        for a in c: inc[a].append(x[i])
    for a in range(V): m.AddExactlyOne(inc[a])
    m.Minimize(sum(x)); s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = workers
    st = s.Solve(m)
    chosen = [cl[i] for i in range(len(cl)) if s.Value(x[i])]
    return chosen, int(s.BestObjectiveBound()), s.StatusName(st)


def alpha(free, tl=60, workers=2):
    """max independent set of the free graph (pairwise-required set): mr >= alpha. returns (value, set, ub, status)"""
    from ortools.sat.python import cp_model
    V = free.shape[0]
    m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in range(V)]
    I, J = np.nonzero(np.triu(free, 1))
    for a, b in zip(I.tolist(), J.tolist()): m.AddBoolOr([x[a].Not(), x[b].Not()])
    m.Maximize(sum(x)); s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = workers
    st = s.Solve(m)
    sel = [a for a in range(V) if s.Value(x[a])]
    return len(sel), sel, int(s.BestObjectiveBound()), s.StatusName(st)


def check_masked(N, free, tol=0):
    """exact/num check that N is a masked identity: diag nonzero, zero on required off-diagonal."""
    V = free.shape[0]
    req = ~free; np.fill_diagonal(req, False)
    off = np.abs(N[req]).max() if req.any() else 0
    dmin = np.abs(np.diag(N)).min()
    return off <= tol and dmin > tol, off, dmin


def rank_modp(M, p=2**31 - 1):
    """exact rank of an integer matrix mod p (gaussian elimination, int64 safe for p < 2^31)."""
    A = np.array(M, dtype=np.int64) % p
    r = 0; rows, cols = A.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c]: piv = i; break
        if piv is None: continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        nz = np.nonzero(A[:, c])[0]
        for i in nz:
            if i != r: A[i] = (A[i] - A[i, c] * A[r]) % p
        r += 1
        if r == rows: break
    return r
