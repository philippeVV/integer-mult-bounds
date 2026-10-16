"""Exact-ish undirected max total multiflow F_* for Tanner sessions (stronger-rank-screens Sec.2).
Sessions = an information set Rs of the code (edges); session k demands from one end of edge e_k to the other through the
internal graph Gamma - Rs (undirected, unit capacities), rate <= 1.  F_* = max total flow.  s/(Wm) >= F_*/W.
Edge-based multicommodity LP (HiGHS).  Certificates: volume bound (E-S)/dist_min (exact) and the LP dual value.
usage: as module: fstar(EV_ends (E x 2), Rs) -> dict"""
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
import networkx as nx

def fstar(ends, Rs, n):
    E = len(ends); Rs = sorted(set(int(x) for x in Rs)); S = len(Rs); isR = np.zeros(E, bool); isR[Rs] = True
    I = np.flatnonzero(~isR); nI = len(I)
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from([tuple(ends[e]) for e in I])
    dists = [nx.shortest_path_length(G, int(ends[e][0]), int(ends[e][1])) if nx.has_path(G, int(ends[e][0]), int(ends[e][1])) else 10**9 for e in Rs]
    dmin = min(dists); vol = nI / dmin
    # variables: for commodity k, arc flows on 2*nI arcs (both directions) + rate t_k
    A = np.concatenate([ends[I], ends[I][:, ::-1]])     # 2nI arcs
    na = len(A); nv = S * na + S
    rows, cols, vals = [], [], []; b_eq = []; r = 0
    for k, e in enumerate(Rs):
        s_, t_ = int(ends[e][0]), int(ends[e][1])
        # conservation at each node: out - in = t_k*(1[s]-1[t])
        for a in range(na):
            u, w = A[a]
            rows += [r + u, r + w]; cols += [k * na + a] * 2; vals += [1.0, -1.0]
        rows += [r + s_, r + t_]; cols += [S * na + k] * 2; vals += [-1.0, 1.0]
        b_eq += [0.0] * n; r += n
    Aeq = sp.csr_matrix((vals, (rows, cols)), shape=(r, nv))
    # capacity: sum_k f_k(a) + f_k(a') <= 1 per undirected edge
    rows, cols = [], []
    for k in range(S):
        for i in range(nI):
            rows += [i, i]; cols += [k * na + i, k * na + nI + i]
    Aub = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(nI, nv))
    c = np.zeros(nv); c[S * na:] = -1.0
    bounds = [(0, None)] * (S * na) + [(0, 1)] * S
    res = linprog(c, A_ub=Aub, b_ub=np.ones(nI), A_eq=Aeq, b_eq=np.array(b_eq), bounds=bounds, method='highs')
    F = -res.fun
    return dict(S=S, E=E, internal=nI, dist_min=int(dmin), dist_mean=float(np.mean(dists)), volume_bound=float(vol),
                F_star=float(F), F_over_W=float(F / S), status=res.status)
