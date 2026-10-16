"""expXH lib: valid-set evaluation for Gram decoders W/2 under parity labels (expNZ/expXA model).
Ports: integer vectors x (rows of X) with form diag(eps), W = X D X^T, W_aa = 2, off-diagonal in {0,+-1}.
Labels B = (W+J) mod 2; h = rank_F2 B; decoder N = W/2 = sum_c eps_c x_c x_c^T / 2 (coordinate terms);
cost C = sum_c rank_F2 B[S_c,S_c], S_c = {a : x_ac != 0}; gamma = C/V."""
import numpy as np


def rank_f2_rows(rows):
    piv = {}
    r = 0
    for v in rows:
        while v:
            lb = v.bit_length() - 1
            if lb in piv:
                v ^= piv[lb]
            else:
                piv[lb] = v; r += 1; break
    return r


def bits_of(M01):
    """rows of a 0/1 numpy matrix -> python ints"""
    out = []
    for row in M01:
        out.append(int(''.join('1' if b else '0' for b in row[::-1]), 2) if len(row) else 0)
    return out


def evaluate(X, eps, check=True):
    X = np.asarray(X, dtype=np.int64)
    V = len(X)
    W = (X * eps) @ X.T
    off = W[~np.eye(V, dtype=bool)]
    if check:
        assert np.all(np.diag(W) == 2), 'norm'
        assert np.all(np.isin(off, (-1, 0, 1))), 'valuation/validity'
    B = ((W + 1) % 2).astype(np.uint8)
    np.fill_diagonal(B, 1)
    h = rank_f2_rows(bits_of(B))
    rk = int(np.linalg.matrix_rank(W.astype(float)))
    C = 0
    det = []
    for c in range(X.shape[1]):
        S = np.nonzero(X[:, c])[0]
        if len(S) == 0:
            continue
        r = rank_f2_rows(bits_of(B[np.ix_(S, S)]))
        det.append(r); C += r
    return dict(V=V, h=h, rank=rk, C=C, gamma=C / V, merit=C / V * (h - 1), det=det)


def max_clique(adj, seconds=60, workers=1, hint=None):
    """CP-SAT maximum clique on a boolean adjacency matrix. returns (clique, optimal?)."""
    from ortools.sat.python import cp_model
    V = len(adj)
    m = cp_model.CpModel()
    x = [m.NewBoolVar('') for _ in range(V)]
    iu = np.argwhere(np.triu(~adj, 1))
    for a, b in iu:
        m.AddBoolOr([x[a].Not(), x[b].Not()])
    m.Maximize(sum(x))
    if hint is not None:
        hs = set(hint)
        for a in range(V):
            m.AddHint(x[a], a in hs)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = seconds; s.parameters.num_workers = workers
    st = s.Solve(m)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return [], False, None
    cl = [a for a in range(V) if s.Value(x[a])]
    cla = np.array(cl, dtype=int)
    if len(cl) and not np.all(adj[np.ix_(cla, cla)] | np.eye(len(cl), dtype=bool)):
        return [], False, None
    return cl, st == cp_model.OPTIMAL, s.BestObjectiveBound()


def clique_ls(adj, iters=20000, restarts=5, seed=0, init=None):
    """Simple tabu/plateau local search for maximum clique (1-swaps). adj: bool numpy matrix (no self loops)."""
    import random
    rng = random.Random(seed)
    V = len(adj)
    A = adj.astype(np.int32)
    best = list(init) if init else []
    for r in range(restarts):
        K = set(best) if (r == 0 and best) else set()
        miss = np.zeros(V, dtype=np.int32)       # number of clique members NOT adjacent (excluding self)
        inK = np.zeros(V, dtype=bool)
        for a in K:
            inK[a] = True
        if K:
            idx = np.array(sorted(K))
            miss = (len(K) - A[:, idx].sum(1)).astype(np.int32)
            miss[idx] -= 0
            for a in idx: miss[a] = 0  # members are adjacent to all others; self not counted (A[a,a]=0 -> miss=1) fix below
            miss[idx] = 0
        tabu = {}
        cur_best = len(K)
        for it in range(iters):
            free = np.nonzero((miss == 0) & ~inK)[0]
            if len(free):
                a = int(free[rng.randrange(len(free))])
                K.add(a); inK[a] = True; miss += (1 - A[:, a]); miss[a] -= 1  # a not adjacent to itself
                miss[a] = 0
            else:
                one = [int(b) for b in np.nonzero((miss == 1) & ~inK)[0] if tabu.get(int(b), -1) < it]
                if not one:
                    if not K: break
                    d = rng.choice(list(K))
                    K.remove(d); inK[d] = False; miss -= (1 - A[:, d]); miss[d] = 0; tabu[d] = it + 7
                    continue
                b = rng.choice(one)
                # remove the unique non-neighbour of b in K
                Kl = np.array(sorted(K)); d = int(Kl[A[b, Kl] == 0][0])
                K.remove(d); inK[d] = False; miss -= (1 - A[:, d]); miss[d] = 0
                tabu[d] = it + 7
                K.add(b); inK[b] = True; miss += (1 - A[:, b]); miss[b] = 0
                miss[np.array(sorted(K))] = 0
            if len(K) > len(best):
                best = sorted(K)
    # verify
    bl = np.array(best)
    assert np.all(adj[np.ix_(bl, bl)] | np.eye(len(bl), dtype=bool))
    return best
