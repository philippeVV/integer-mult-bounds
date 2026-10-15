"""Certificate-emitting zero-robust propagation engine (see engine.py docstring).

propagate(normals, caps) -> ("SURVIVE", state, steps) or ("INFEASIBLE", final, steps)
steps: list of kill records, each with its justification:
  ("L2", g, side, r, Q)            side in {"R","C"}; Q = tuple of live positions (|Q| = cap_g)
  ("L3", g, side, r, x, y, walk)   walk = odd-length walk x..y in the owned graph of g (cap_g == 2)
final contradiction:
  ("L0", r) | ("L1", Q) | ("L2rank", g, Q) | ("L3odd", g, closed_walk)
All positions are live positions 0..111 (source index = LIVE[pos]).
"""
import numpy as np
from model import NL, ONB, COSTAR, popc, bits, LIVE

_CL = None
def load_cliques():
    global _CL
    if _CL is None:
        import pickle
        cl = pickle.load(open("cliques.pkl", "rb"))
        M = np.zeros((len(cl), NL), dtype=np.float32)
        for i, q in enumerate(cl):
            for b in bits(q):
                M[i, b] = 1
        _CL = (cl, M, M.sum(axis=1).astype(np.int32))
    return _CL

def vecs(masks):
    V = np.zeros((NL, len(masks)), dtype=np.float32)
    for j, m in enumerate(masks):
        for b in bits(m):
            V[b, j] = 1
    return V

def owners_tables(R, C):
    rowg = [0] * NL
    colg = [0] * NL
    for g in range(len(R)):
        for b in bits(R[g]):
            rowg[b] |= 1 << g
        for b in bits(C[g]):
            colg[b] |= 1 << g
    return rowg, colg

def l1_check(R, C, caps):
    cl, M, sz = load_cliques()
    G = len(R)
    A = M @ vecs(R + C)          # cliques x 2G
    A = A.astype(np.int32)
    tot = np.zeros(len(cl), dtype=np.int32)
    for g in range(G):
        tot += np.minimum(np.minimum(A[:, g], A[:, G + g]), caps[g])
    bad = np.nonzero(tot < sz)[0]
    if len(bad):
        i = int(bad[np.argmin(tot[bad] - sz[bad])])
        return tuple(bits(cl[i]))
    return None

def owned_graph(g, R, C, rowg, colg):
    me = 1 << g
    def alone(r, s):
        return (rowg[r] & colg[s]) & ~me == 0
    act = R[g] & C[g]
    V = [x for x in bits(act) if alone(x, x)]
    Vm = 0
    for x in V:
        Vm |= 1 << x
    adj = {}
    for x in V:
        m = 0
        for y in bits(ONB[x] & Vm):
            if alone(x, y) and alone(y, x):
                m |= 1 << y
        adj[x] = m
    return V, Vm, adj, alone

def propagate(normals, caps, maxiter=500, use_l3=True):
    G = len(normals)
    R = [COSTAR[a] for a in normals]
    C = [COSTAR[a] for a in normals]
    caps = list(caps)
    steps = []
    for it in range(maxiter):
        rowg, colg = owners_tables(R, C)
        for r in range(NL):
            if rowg[r] & colg[r] == 0:
                return "INFEASIBLE", ("L0", r), steps
        Q = l1_check(R, C, caps)
        if Q is not None:
            return "INFEASIBLE", ("L1", Q), steps
        newR = list(R)
        newC = list(C)
        new_steps = []
        for g in range(G):
            c = caps[g]
            V, Vm, adj, alone = owned_graph(g, R, C, rowg, colg)
            # L2: cliques of size c in owned graph
            Qs = []
            def ext(cur, cand, size):
                if size == c:
                    Qs.append(cur)
                    return
                cc = cand
                while cc:
                    bb = cc & -cc
                    i = bb.bit_length() - 1
                    cc ^= bb
                    ext(cur | bb, cc & adj[i], size + 1)
            for x in V:
                ext(1 << x, adj[x] & ~((1 << (x + 1)) - 1), 1)
            for Qm in Qs:
                qs = bits(Qm)
                cand = Vm
                for q in qs:
                    cand &= adj[q]
                if cand:
                    extra = (cand & -cand).bit_length() - 1
                    return "INFEASIBLE", ("L2rank", g, tuple(sorted(qs + [extra]))), steps + new_steps
                common = (1 << NL) - 1
                for q in qs:
                    common &= ONB[q]
                for r in bits(common & newR[g]):
                    if all(alone(r, q) for q in qs):
                        newR[g] &= ~(1 << r)
                        new_steps.append(("L2", g, "R", r, tuple(qs)))
                for r in bits(common & newC[g]):
                    if all(alone(q, r) for q in qs):
                        newC[g] &= ~(1 << r)
                        new_steps.append(("L2", g, "C", r, tuple(qs)))
            # L3 for capacity 2
            if use_l3 and c == 2:
                color = {}
                parent = {}
                for x0 in V:
                    if x0 in color or adj[x0] == 0:
                        continue
                    color[x0] = 0
                    parent[x0] = None
                    comp = [x0]
                    st = [x0]
                    while st:
                        u = st.pop(0)
                        for w in bits(adj[u]):
                            if w not in color:
                                color[w] = 1 - color[u]
                                parent[w] = u
                                comp.append(w)
                                st.append(w)
                            elif color[w] == color[u]:
                                def up(z):
                                    p = [z]
                                    while parent[p[-1]] is not None:
                                        p.append(parent[p[-1]])
                                    return p
                                walk = list(reversed(up(u))) + up(w)   # root..u, w..root : odd closed walk
                                return "INFEASIBLE", ("L3odd", g, tuple(walk)), steps + new_steps
                    def walk_between(x, y):
                        def up(z):
                            p = [z]
                            while parent[p[-1]] is not None:
                                p.append(parent[p[-1]])
                            return p
                        return tuple(list(up(x)) + list(reversed(up(y)))[1:])   # x..root..y
                    Xs = [u for u in comp if color[u] == 0]
                    Ys = [u for u in comp if color[u] == 1]
                    Xm = sum(1 << u for u in Xs)
                    Ym = sum(1 << u for u in Ys)
                    for r in bits(newR[g]):
                        xs = [x for x in bits(ONB[r] & Xm) if alone(r, x)]
                        ys = [y for y in bits(ONB[r] & Ym) if alone(r, y)]
                        if xs and ys:
                            newR[g] &= ~(1 << r)
                            new_steps.append(("L3", g, "R", r, xs[0], ys[0], walk_between(xs[0], ys[0])))
                    for r in bits(newC[g]):
                        xs = [x for x in bits(ONB[r] & Xm) if alone(x, r)]
                        ys = [y for y in bits(ONB[r] & Ym) if alone(y, r)]
                        if xs and ys:
                            newC[g] &= ~(1 << r)
                            new_steps.append(("L3", g, "C", r, xs[0], ys[0], walk_between(xs[0], ys[0])))
        if not new_steps:
            return "SURVIVE", (R, C), steps
        steps += new_steps
        R, C = newR, newC
    return "SURVIVE", (R, C), steps
