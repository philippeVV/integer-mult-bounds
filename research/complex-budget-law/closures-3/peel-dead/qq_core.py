"""expQQ core: XOR (GL(2,F2)) circuits, physical graph, union-find with integer offsets,
frame problems restricted to <=2 frame classes per role route.

Conventions follow upstream scripts/finite_bit_contract.py:
  nodes: source w -> gates (time order) -> sink w on every role route; physical edges are consecutive pairs;
  endpoint identity  T_{rho(w)} - S_w = I  (T = sink frame, S = source frame);
  s = sum over physical edges of rank(frame(head) - frame(tail)).
Frames that are equal up to a multiple of I are tracked by union-find with integer offsets (gain graph over Z).
"""
from itertools import combinations, product

# GL(2,F2) acting on the pair of values (va, vb) of roles (a, b)
OPS = (
    lambda a, b: (a, b),          # identity (still forces a common frame on both ports)
    lambda a, b: (b, a),          # swap
    lambda a, b: (a ^ b, b),      # a += b
    lambda a, b: (a, b ^ a),      # b += a
    lambda a, b: (b, a ^ b),      # 3-cycle
    lambda a, b: (a ^ b, a),      # 3-cycle inverse
)
OP_XORS = ([], [[0, 1], [1, 0], [0, 1]], [[0, 1]], [[1, 0]], None, None)


def run(W, pairs, ops):
    v = [1 << w for w in range(W)]
    for (a, b), o in zip(pairs, ops):
        v[a], v[b] = OPS[o](v[a], v[b])
    return v


def perm_of(W, vals):
    """rho(input)=output, or None if the F2 map is not a permutation matrix."""
    if sorted(vals) != [1 << w for w in range(W)]:
        return None
    rho = [None] * W
    for out, val in enumerate(vals):
        rho[val.bit_length() - 1] = out
    return tuple(rho)


def achievable_perms(W, pairs):
    """All rho realised by some choice of the 6 ops per gate (BFS over distinct value tuples)."""
    states = {tuple(1 << w for w in range(W))}
    for a, b in pairs:
        nxt = set()
        for st in states:
            for f in OPS:
                v = list(st)
                v[a], v[b] = f(v[a], v[b])
                nxt.add(tuple(v))
        states = nxt
    out = set()
    for st in states:
        r = perm_of(W, st)
        if r is not None:
            out.add(r)
    return out


def routes(W, pairs):
    """node ids: sources 0..W-1, gates W..W+n-1, sinks W+n..W+n+W-1. Returns per-role node list."""
    n = len(pairs)
    R = [[w] for w in range(W)]
    for i, (a, b) in enumerate(pairs):
        R[a].append(W + i)
        R[b].append(W + i)
    for w in range(W):
        R[w].append(W + n + w)
    return R


def routable(W, pairs, rho):
    """Token routing: each gate either passes or swaps its two tokens."""
    want = [None] * W
    for w, o in enumerate(rho):
        want[o] = w
    states = {tuple(range(W))}
    for a, b in pairs:
        nxt = set()
        for st in states:
            nxt.add(st)
            l = list(st); l[a], l[b] = l[b], l[a]
            nxt.add(tuple(l))
        states = nxt
    return tuple(want) in states


class UF:
    """union-find with integer offsets: value(x) = value(root(x)) + off(x)."""
    __slots__ = ('p', 'o')

    def __init__(self, N):
        self.p = list(range(N)); self.o = [0] * N

    def find(self, x):
        p, o = self.p, self.o
        path = []
        while p[x] != x:
            path.append(x); x = p[x]
        # compress
        acc = 0
        for y in reversed(path):
            acc += o[y]
            o[y] = acc
            p[y] = x
        return x

    def off(self, x):
        self.find(x)
        return self.o[x] if self.p[x] != x else 0

    def union(self, x, y, g):
        """impose value(y) = value(x) + g. Return False on contradiction."""
        rx, ry = self.find(x), self.find(y)
        ox = self.o[x] if x != rx else 0
        oy = self.o[y] if y != ry else 0
        if rx == ry:
            return oy - ox == g
        # value(ry) = value(y) - oy = value(x) + g - oy = value(rx) + ox + g - oy
        self.p[ry] = rx; self.o[ry] = ox + g - oy
        return True


def scalar_min_changes(W, pairs, rho, limit=None):
    """Scalar (m=1) frames, any number of climbs: minimum number of physical edges with a nonzero change.
    Returns the minimum if < limit (default W), else returns limit (meaning >= limit). Exact."""
    if limit is None:
        limit = W
    R = routes(W, pairs)
    n = len(pairs)
    N = 2 * W + n
    edges = [(r[i], r[i + 1]) for r in R for i in range(len(r) - 1)]
    for k in range(0, limit):
        for drop in combinations(range(len(edges)), k):
            ds = set(drop)
            uf = UF(N)
            ok = True
            for w in range(W):
                if not uf.union(w, W + n + rho[w], 1):
                    ok = False; break
            if not ok:
                continue
            for j, (u, v) in enumerate(edges):
                if j in ds:
                    continue
                if not uf.union(u, v, 0):
                    ok = False; break
            if ok:
                return k
    return limit


def route_patterns(k, monotone):
    """class sequences (0=A, 1=B) for the k+2 nodes of a route with k gates; source is class A."""
    L = k + 2
    if monotone:
        for j in range(L + 1):          # nodes with index >= j are B  (j = L: no climb)
            if j == 0:
                continue
            yield tuple(0 if i < j else 1 for i in range(L))
    else:
        for bits in product((0, 1), repeat=L - 1):
            yield (0,) + bits


def solve_pattern(W, n, R, rho, pat):
    """Exact gluing for one two-class pattern. Returns None if the hard equalities are inconsistent, else
    (hard_cost_units, soft) where soft = list of (rootA, rootB, offset, weight): the class difference
    B_w - A_w equals P_rootB - P_rootA + offset*I and is paid weight times."""
    N = 2 * W + n
    uf = UF(N + 2 * W)
    for w in range(W):
        if not uf.union(w, W + n + rho[w], 1):
            return None
    for w in range(W):
        for node, c in zip(R[w], pat[w]):
            if not uf.union(N + 2 * w + c, node, 0):
                return None
    hard = 0
    soft = []
    for w in range(W):
        p = pat[w]
        c = sum(p[i] != p[i + 1] for i in range(len(p) - 1))
        if c == 0:
            continue
        a, b = N + 2 * w, N + 2 * w + 1
        ra, rb = uf.find(a), uf.find(b)
        oa, ob = uf.off(a), uf.off(b)
        if ra == rb:
            if ob - oa != 0:
                hard += c
        else:
            soft.append((ra, rb, ob - oa, c))
    return hard, soft


def soft_scalar_min(soft):
    """min total weight of soft relations left unsatisfied (scalar potentials, integer gains)."""
    best = sum(s[3] for s in soft)
    k = len(soft)
    for mask in range(1 << k):
        cost = sum(soft[i][3] for i in range(k) if not mask >> i & 1)
        if cost >= best:
            continue
        uf = {}
        def find(x):
            # tiny dict-based offset union-find
            path = []
            while x in uf and uf[x][0] != x:
                path.append(x); x = uf[x][0]
            return x
        # simple: rebuild with UF over relabelled roots
        roots = sorted({s[0] for s in soft} | {s[1] for s in soft})
        idx = {r: i for i, r in enumerate(roots)}
        u = UF(len(roots)); ok = True
        for i in range(k):
            if mask >> i & 1:
                ra, rb, g, _ = soft[i]
                if not u.union(idx[ra], idx[rb], g):
                    ok = False; break
        if ok:
            best = cost
    return best


def soft_cycles(soft):
    """all simple unbalanced cycles (as index tuples) of the soft multigraph (no loops: roots differ)."""
    k = len(soft)
    out = []
    for mask in range(1, 1 << k):
        idxs = [i for i in range(k) if mask >> i & 1]
        deg = {}
        for i in idxs:
            for v in soft[i][:2]:
                deg[v] = deg.get(v, 0) + 1
        if any(d != 2 for d in deg.values()):
            continue
        # walk the cycle to check connectivity and get the signed gain
        start = soft[idxs[0]][0]
        used = set(); cur = start; gain = 0
        while True:
            nxt = None
            for i in idxs:
                if i in used:
                    continue
                a, b, g, _ = soft[i]
                if a == cur:
                    nxt = (i, b, g); break
                if b == cur:
                    nxt = (i, a, -g); break
            if nxt is None:
                break
            used.add(nxt[0]); cur = nxt[1]; gain += nxt[2]
        if len(used) == len(idxs) and cur == start and gain != 0:
            out.append(tuple(idxs))
    return out


def soft_lp(soft):
    """fractional packing of unbalanced cycles with capacity = weight (lower bound on matrix soft cost / m)."""
    cyc = soft_cycles(soft)
    if not cyc:
        return 0.0
    from scipy.optimize import linprog
    import numpy as np
    A = np.zeros((len(soft), len(cyc)))
    for j, C in enumerate(cyc):
        for i in C:
            A[i, j] = 1
    res = linprog(-np.ones(len(cyc)), A_ub=A, b_ub=[s[3] for s in soft], bounds=(0, None), method='highs')
    return -res.fun
