"""Residual routing: route every unmatched demand unit (a->b) along an undirected path of residual edges
(unmatched climbs + E-arcs), edge-disjoint (each physical edge capacity 1). Integral backtracking.
A success + direct routing of matched pairs = W edge-disjoint own-pair paths in the class graph, hence
s >= W m for every frame assignment (stronger-rank-screens §2, applied to the class graph)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir


def residual_routable(levels, arcs, d, c):
    edges = []
    for (a, b), m in c.items():
        edges += [(a, b)] * m
    edges += list(arcs)
    dem = []
    for (a, b), m in d.items():
        dem += [(a, b)] * m
    adj = {}
    for i, (x, y) in enumerate(edges):
        adj.setdefault(x, []).append((i, y)); adj.setdefault(y, []).append((i, x))

    def paths(s, t, used):
        # simple paths s->t avoiding used edges
        out = []
        stack = [(s, [], {s})]
        while stack:
            v, pe, seen = stack.pop()
            if v == t:
                out.append(pe); continue
            for i, w in adj.get(v, []):
                if i in used or i in pe or w in seen:
                    continue
                stack.append((w, pe + [i], seen | {w}))
        return out

    def bt(j, used):
        if j == len(dem):
            return True
        s, t = dem[j]
        for p in paths(s, t, used):
            if bt(j + 1, used | set(p)):
                return True
        return False
    return bt(0, frozenset())
