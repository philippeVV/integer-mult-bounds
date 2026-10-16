"""Baer-quotient layouts for PG(2,q)-M_{d0}, q=4^f.  Singer model: point i ~ line i+d, d in D' = D-{d0}.
Baer classes = cosets of the order-beta subgroup (beta=q+sqrt q+1): class(i) = i mod c, c = q-sqrt q+1.
D mod c = {r0 x (sqrt q+1), every other residue x1}.  Removing d0 with r1 = d0-r0 != 0 mod c kills exactly the
class pairs (a, a+r1) -- so arranging classes so that a+r1 lies 'across' every prefix makes all M-edges cut edges.
Layout = class sequence (t -> x0 + s*t mod c) x within-class order.  Exact |R| via expNI mcf_eval (K = q/2)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, math, numpy as np
from ortools.graph.python import min_cost_flow

def singer_D(e):
    from dset import dset
    D = dset(e); q = 2 ** e; return q * q + q + 1, D

class G:
    def __init__(s, q, n, D, d0):
        s.q, s.n, s.D, s.d0 = q, n, list(D), d0
        s.r = math.isqrt(q); s.c = q - s.r + 1; s.beta = q + s.r + 1; assert s.c * s.beta == n
        from collections import Counter
        cnt = Counter(d % s.c for d in D); s.r0 = max(cnt, key=cnt.get); assert cnt[s.r0] == s.r + 1
        s.r1 = (d0 - s.r0) % s.c
        s.Dp = [d for d in D if d != d0]; s.V = 2 * n
        s.edges = np.array([(i, n + (i + d) % n) for i in range(n) for d in s.Dp], dtype=np.int64)
        nb = np.zeros((s.V, q), dtype=np.int64); fill = np.zeros(s.V, dtype=np.int64)
        for a, b in s.edges:
            nb[a, fill[a]] = b; fill[a] += 1; nb[b, fill[b]] = a; fill[b] += 1
        assert (fill == q).all(); s.nb = nb
    def cls(s, v):  # Baer class of a vertex: points i -> i mod c, lines j -> (j - r0) mod c (Baer lines of class x)
        return v % s.c if v < s.n else (v - s.n - s.r0) % s.c
    def profile(s, order):
        pos = np.empty(s.V, dtype=np.int64); pos[np.asarray(order)] = np.arange(s.V)
        a, b = s.edges[:, 0], s.edges[:, 1]
        lo = np.minimum(pos[a], pos[b]); hi = np.maximum(pos[a], pos[b])
        diff = np.zeros(s.V + 1, dtype=np.int64); np.add.at(diff, lo, 1); np.add.at(diff, hi, -1)
        return np.cumsum(diff)[:-1], pos   # cut after position t
    def mcf(s, pos, want=False):
        K = s.q // 2; BIG = 10 ** 6; V = s.V
        m = min_cost_flow.SimpleMinCostFlow()
        a, b = s.edges[:, 0], s.edges[:, 1]
        fw = pos[a] < pos[b]; u = np.where(fw, a, b); w = np.where(fw, b, a)
        arcs = m.add_arcs_with_capacity_and_unit_cost(u, w, np.ones(len(u), np.int64), np.ones(len(u), np.int64))
        indeg = np.bincount(w, minlength=V); sup = K - indeg
        for v in range(V):
            m.set_node_supply(v, int(sup[v]))
            if sup[v] > 0: m.add_arc_with_capacity_and_unit_cost(v, V, int(sup[v]), BIG)
            elif sup[v] < 0: m.add_arc_with_capacity_and_unit_cost(V, v, int(-sup[v]), BIG)
        assert m.solve() == m.OPTIMAL
        cost = m.optimal_cost(); slack = cost // (2 * BIG); R = cost - slack * 2 * BIG
        if want:
            fl = np.array([m.flow(x) for x in arcs]); return slack, R, fl, u, w
        return slack, R

def class_seq(g, s_, x0=0):
    return [(x0 + s_ * t) % g.c for t in range(g.c)]

def members(g, x):
    pts = [x + g.c * t for t in range(g.beta)]
    lns = [g.n + (x + g.r0 + g.c * t) % g.n for t in range(g.beta)]
    return pts, lns
