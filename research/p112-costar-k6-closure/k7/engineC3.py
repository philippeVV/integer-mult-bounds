"""engineC: engine3 with the clique scan (L1 / L2G tightness) done in C (cscan.so); same rules and step format."""

import numpy as np, ctypes, os, pickle
_HERE = os.path.dirname(os.path.abspath(__file__))
_LIB = ctypes.CDLL(os.path.join(_HERE, "cscan.so"))
_MAXT = 400000
_TIDX = (ctypes.c_int * _MAXT)(); _TGAM = (ctypes.c_uint32 * _MAXT)()
_CD = None
def _cdata():
    global _CD
    if _CD is None:
        cl = pickle.load(open(os.path.join(_HERE, "cliques.pkl"), "rb"))
        M64 = (1 << 64) - 1
        lo = np.array([q & M64 for q in cl], dtype=np.uint64); hi = np.array([q >> 64 for q in cl], dtype=np.uint64)
        sz = np.array([q.bit_count() for q in cl], dtype=np.int8)
        _CD = (cl, lo.ctypes.data_as(ctypes.c_void_p), hi.ctypes.data_as(ctypes.c_void_p), sz.ctypes.data_as(ctypes.c_void_p), lo, hi, sz)
    return _CD[:4]
from model import NL, ONB, COSTAR, popc, bits, LIVE
from engine2 import load_cliques, vecs, owners_tables, owned_graph

def propagate(normals, caps, maxiter=500, use_l3=True, R0=None, C0=None):
    cl, QLO, QHI, QSZ = _cdata()
    G = len(normals)
    R = list(R0) if R0 is not None else [COSTAR[a] for a in normals]
    C = list(C0) if C0 is not None else [COSTAR[a] for a in normals]
    caps = list(caps)
    capc = (ctypes.c_int * G)(*caps)
    steps = []
    M64 = (1 << 64) - 1
    for it in range(maxiter):
        rowg, colg = owners_tables(R, C)
        for r in range(NL):
            if rowg[r] & colg[r] == 0:
                return "INFEASIBLE", ("L0", r), steps
        Rlo = (ctypes.c_uint64 * G)(*[x & M64 for x in R]); Rhi = (ctypes.c_uint64 * G)(*[x >> 64 for x in R])
        Clo = (ctypes.c_uint64 * G)(*[x & M64 for x in C]); Chi = (ctypes.c_uint64 * G)(*[x >> 64 for x in C])
        viol = ctypes.c_int(-1)
        nt = _LIB.scan(len(cl), QLO, QHI, QSZ, G, Rlo, Rhi, Clo, Chi, capc, ctypes.byref(viol), _TIDX, _TGAM, _MAXT)
        if viol.value >= 0:
            return "INFEASIBLE", ("L1", tuple(bits(cl[viol.value]))), steps
        if nt > _MAXT:
            raise RuntimeError("too many tight cliques")
        newR = list(R)
        newC = list(C)
        new_steps = []
        # ---- L2G ----
        for j in range(nt):
            Qm = cl[_TIDX[j]]
            qs = bits(Qm)
            gm = _TGAM[j]
            gam = [g for g in range(G) if gm >> g & 1]
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
