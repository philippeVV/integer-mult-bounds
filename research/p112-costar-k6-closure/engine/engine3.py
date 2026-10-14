"""engine2 + L2G (tight multi-group clique kill), vectorized over all cliques.

L2G: Q orthogonal clique; Gamma(Q) = groups alive on some row of Q AND some column of Q
(exactly the groups that can contribute to the Q x Q block).  If sum_{a in Gamma} c_a = |Q|
then D_Gamma[Q], T_Gamma[Q] are square and invertible (their product is 4I).  For r notin Q
orthogonal to all of Q such that no group outside Gamma is alive at (r, q) for any q
(row r in R_b and q in C_b), D_a[r] = 0 for all a in Gamma.  Symmetric for columns.
(L2 = the case Gamma = {a}.)
Certificate step: ("L2G", gamma_tuple, side, r, Q)
"""
import numpy as np
from model import NL, ONB, COSTAR, popc, bits, LIVE
from engine2 import load_cliques, vecs, owners_tables, owned_graph

def propagate(normals, caps, maxiter=500, use_l3=True, R0=None, C0=None):
    cl, M, sz = load_cliques()
    G = len(normals)
    R = list(R0) if R0 is not None else [COSTAR[a] for a in normals]
    C = list(C0) if C0 is not None else [COSTAR[a] for a in normals]
    caps = list(caps)
    capv = np.array(caps, dtype=np.int32)
    steps = []
    for it in range(maxiter):
        rowg, colg = owners_tables(R, C)
        for r in range(NL):
            if rowg[r] & colg[r] == 0:
                return "INFEASIBLE", ("L0", r), steps
        A = (M @ vecs(R + C)).astype(np.int32)
        AR, AC = A[:, :G], A[:, G:]
        b = np.minimum(np.minimum(AR, AC), capv)
        tot = b.sum(axis=1)
        bad = np.nonzero(tot < sz)[0]
        if len(bad):
            i = int(bad[np.argmin(tot[bad] - sz[bad])])
            return "INFEASIBLE", ("L1", tuple(bits(cl[i]))), steps
        newR = list(R)
        newC = list(C)
        new_steps = []
        # ---- L2G ----
        touch = (AR > 0) & (AC > 0)                    # cliques x groups
        ctouch = (touch * capv).sum(axis=1)
        tight = np.nonzero(ctouch == sz)[0]
        for i in tight:
            Qm = cl[int(i)]
            qs = bits(Qm)
            gam = [g for g in range(G) if touch[i, g]]
            gm = sum(1 << g for g in gam)
            common = (1 << NL) - 1
            for q in qs:
                common &= ONB[q]
            colQ = 0      # groups alive on some column of Q
            rowQ = 0
            for q in qs:
                colQ |= colg[q]
                rowQ |= rowg[q]
            outside_col = colQ & ~gm
            outside_row = rowQ & ~gm
            for r in bits(common):
                # row kill: no group outside gamma with r in R and some q in C
                if rowg[r] & outside_col == 0:
                    hit = [g for g in gam if newR[g] >> r & 1]
                    if hit:
                        for g in hit:
                            newR[g] &= ~(1 << r)
                        new_steps.append(("L2G", tuple(gam), "R", r, tuple(qs)))
                if colg[r] & outside_row == 0:
                    hit = [g for g in gam if newC[g] >> r & 1]
                    if hit:
                        for g in hit:
                            newC[g] &= ~(1 << r)
                        new_steps.append(("L2G", tuple(gam), "C", r, tuple(qs)))
        # ---- L3 (capacity-2 bipartite lemma) ----
        if use_l3:
            for g in range(G):
                if caps[g] != 2:
                    continue
                V, Vm, adj, alone = owned_graph(g, R, C, rowg, colg)
                color = {}
                parent = {}
                def up(z):
                    p = [z]
                    while parent[p[-1]] is not None:
                        p.append(parent[p[-1]])
                    return p
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
                                walk = list(reversed(up(u))) + up(w)
                                return "INFEASIBLE", ("L3odd", g, tuple(walk)), steps + new_steps
                    Xm = sum(1 << u for u in comp if color[u] == 0)
                    Ym = sum(1 << u for u in comp if color[u] == 1)
                    for r in bits(newR[g]):
                        xs = [x for x in bits(ONB[r] & Xm) if alone(r, x)]
                        ys = [y for y in bits(ONB[r] & Ym) if alone(r, y)]
                        if xs and ys:
                            newR[g] &= ~(1 << r)
                            new_steps.append(("L3", g, "R", r, xs[0], ys[0], tuple(up(xs[0]) + list(reversed(up(ys[0])))[1:])))
                    for r in bits(newC[g]):
                        xs = [x for x in bits(ONB[r] & Xm) if alone(x, r)]
                        ys = [y for y in bits(ONB[r] & Ym) if alone(y, r)]
                        if xs and ys:
                            newC[g] &= ~(1 << r)
                            new_steps.append(("L3", g, "C", r, xs[0], ys[0], tuple(up(xs[0]) + list(reversed(up(ys[0])))[1:])))
        if not new_steps:
            return "SURVIVE", (R, C), steps
        steps += new_steps
        R, C = newR, newC
    return "SURVIVE", (R, C), steps
