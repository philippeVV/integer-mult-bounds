"""Exhaustive (solver-free) laminar-order search by DFS over edge prefixes with frame-state memoisation.
Used to certify obstructions independently of CP-SAT.  State = tuple of frame sets of all wires + used edges."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys
sys.setrecursionlimit(10000)

def dfs_laminar(B, edges, x0=None, t0=None):
    srcs = sorted({b for b, a in edges}); tgts = sorted({a for b, a in edges})
    xs0 = {b: frozenset({b} | set((x0 or {}).get(b, ()))) for b in srcs}
    ts0 = {a: frozenset((t0 or {}).get(a, ())) for a in tgts}
    seen = set()
    E = len(edges)
    def rec(used, xs, ts, order):
        if used == (1 << E) - 1:
            return order
        key = (used, tuple(xs[b] for b in srcs), tuple(ts[a] for a in tgts))
        if key in seen: return None
        seen.add(key)
        for i in range(E):
            if used >> i & 1: continue
            b, a = edges[i]
            U = xs[b] | ts[a]
            if any(B[c][a] for c in U): continue
            xs2 = dict(xs); ts2 = dict(ts); xs2[b] = U; ts2[a] = U
            r = rec(used | 1 << i, xs2, ts2, order + [i])
            if r is not None: return r
        return None
    res = rec(0, xs0, ts0, [])
    return res, len(seen)

if __name__ == '__main__':
    # abstract bowtie: centre 0, triangles {0,1,2} and {0,3,4}; free pairs inside triangles, cross pairs required
    n = 5; free = {(0, 1), (0, 2), (1, 2), (0, 3), (0, 4), (3, 4)}
    B = [[1 if i == j or ((min(i, j), max(i, j)) not in free) else 0 for j in range(n)] for i in range(n)]
    E = [(a, b) for a, b in free] + [(b, a) for a, b in free]
    r, st = dfs_laminar(B, E)
    print('bowtie (all 12 directed corrections):', 'LAMINAR' if r else 'NO laminar order', 'states', st)
    # directed sub-obstructions of the bowtie: remove each edge
    crit = []
    for i in range(len(E)):
        E2 = E[:i] + E[i + 1:]
        r2, _ = dfs_laminar(B, E2)
        crit.append(r2 is not None)
    print('removing any single correction makes it laminar:', all(crit), crit)
    # single triangle with all 6 corrections + pendant required ports: laminar?
    tri = [(0, 1), (0, 2), (1, 2)]
    B3 = [[1 if i == j else 0 for j in range(3)] for i in range(3)]
    E3 = [(a, b) for a, b in tri] + [(b, a) for a, b in tri]
    print('free triangle:', 'LAMINAR' if dfs_laminar(B3, E3)[0] else 'NO')
    # one-way bowtie (antisymmetric orientation) for asymmetric decoders
    import itertools
    best = None
    for mask in range(1 << 6):
        Ed = [(a, b) if mask >> k & 1 else (b, a) for k, (a, b) in enumerate(sorted(free))]
        if dfs_laminar(B, Ed)[0] is None:
            best = Ed; break
    print('some one-way orientation of the bowtie non-laminar:', best)
    print('PASS_DFS')
