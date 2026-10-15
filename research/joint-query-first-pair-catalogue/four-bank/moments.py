#!/usr/bin/env python3
"""Exact spectral moments tr(T^k)/N, k<=4, of T = sum_{s in S} B_s* B_s / 48 for stage triples S,
on the actual QQ16 routes (regular cover, payload-constant fields), by explicit group-algebra
block contraction C=(B_s B_t*)_{s,t}: tr C^2 = sum ||C||^2, tr C^3 = <C^2, C>, tr C^4 = ||C^2||^2.
No collision classification is assumed.  Also the rank bound rank T >= tr p(T),
p(t) = 1-(1-alpha t-beta t^2)^2, and the erased-family defects.
usage: moments.py ARCHIVE OUT.json [conv]   conv in {hom, anti}"""
import json, sys, time
from fractions import Fraction as Q
sys.path.insert(0, __file__.rsplit('/', 1)[0] + '/..')
import pair_catalogue as pc

def main():
    t0 = time.time()
    local = pc.load(sys.argv[1])
    cols = pc.query_columns(local)
    ph = local["DATA-PHASES.json"]
    bases = [tuple(r["idle_children"][0]["orthogonal_axes"][:8]) for r in ph["recipes"]]
    conv = sys.argv[3] if len(sys.argv) > 3 else 'hom'
    U = pc.U
    rr = {p: [cols[p][i] * cols[p][j] for i in range(8) for j in range(8)] for p in U}
    res = {}
    import itertools
    trip = [(0, 1, 2), (2, 3, 4), (4, 5, 6), (6, 7, 8)]
    if len(sys.argv) > 4:
        w, nw = map(int, sys.argv[4].split('/'))
        trip = [S for i, S in enumerate(itertools.combinations(range(9), 3)) if i % nw == w]
    for S in trip:
        coords = pc.coords_for(list(S))
        ID = tuple(1 << i for i in range(len(coords)))
        g = {s: [pc.stage_route(bases[p], s, coords) for p in range(112)] for s in S}
        for s in S:
            for p in U:
                pc.need(pc.mul(g[s][p], g[s][p]) == ID, 'involution')
        def rel(s, t, p):  # element of B_s B_t^* at port p (routes are involutions)
            return pc.mul(g[s][p], g[t][p]) if conv == 'hom' else pc.mul(g[t][p], g[s][p])
        Cb = {}
        for s in S:
            for t in S:
                d = {}
                if s == t:
                    d[ID] = [48 * (i == j) for i in range(8) for j in range(8)]
                else:
                    for p in U:
                        x = rel(s, t, p)
                        pc.need(x not in d, 'relative routes distinct')
                        d[x] = list(rr[p])
                Cb[s, t] = d
        def mm(A, B):
            return [sum(A[8 * i + k] * B[8 * k + j] for k in range(8)) for i in range(8) for j in range(8)]
        def prod(X, Y):
            out = {}
            for gx, mx in X.items():
                for gy, my in Y.items():
                    z = pc.mul(gx, gy) if conv == 'hom' else pc.mul(gy, gx)
                    m = mm(mx, my)
                    o = out.get(z)
                    if o is None:
                        out[z] = m
                    else:
                        for i in range(64):
                            o[i] += m[i]
            return {k: v for k, v in out.items() if any(v)}
        C2 = {}
        for s in S:
            for t in S:
                acc = {}
                for u in S:
                    for k, v in prod(Cb[s, u], Cb[u, t]).items():
                        o = acc.setdefault(k, [0] * 64)
                        for i in range(64):
                            o[i] += v[i]
                C2[s, t] = {k: v for k, v in acc.items() if any(v)}
        tr1 = sum(Cb[s, s][ID][9 * i] for s in S for i in range(8))
        tr2 = sum(sum(x * x for x in v) for blk in Cb.values() for v in blk.values())
        tr3 = sum(sum(a * b for a, b in zip(v, Cb[s, t].get(k, [0] * 64))) for (s, t), blk in C2.items() for k, v in blk.items())
        tr4 = sum(sum(x * x for x in v) for blk in C2.values() for v in blk.values())
        # check tr C^2 also equals identity coefficient of C^2 diagonal blocks
        pc.need(tr2 == sum(C2[s, s].get(ID, [0] * 64)[9 * i] for s in S for i in range(8)), 'tr C^2 consistency')
        m = [Q(tr1, 48), Q(tr2, 48 ** 2), Q(tr3, 48 ** 3), Q(tr4, 48 ** 4)]
        det = m[1] * m[3] - m[2] * m[2]
        al = (m[0] * m[3] - m[1] * m[2]) / det
        be = (m[1] * m[1] - m[0] * m[2]) / det
        rank_lb = 2 * al * m[0] + 2 * be * m[1] - (al * al * m[1] + 2 * al * be * m[2] + be * be * m[3])
        # sanity: p(t) <= 1 everywhere is automatic; rank bound via second moment too
        res[str(S)] = dict(moments=[str(x) for x in m], alpha=str(al), beta=str(be), rank_lb_over_N=str(rank_lb),
                           rank_lb_2nd_moment=str(m[0] ** 2 / m[1]),
                           C2_terms=sum(len(b) for b in C2.values()))
        print(S, 'moments', [str(x) for x in m], 'rank>=', rank_lb, float(rank_lb), f'{time.time()-t0:.0f}s', flush=True)
    json.dump(res, open(sys.argv[2], 'w'), indent=1)

main()
