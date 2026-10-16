"""Leaf-extension variants of the BH seed (all leaves capacity 1)."""
import sys; sys.path.insert(0, __file__.rsplit('/',1)[0])
from bh import CORE, SESS
from collections import defaultdict
def full():
    E = list(CORE); nv = 182; ss = []
    for q, s, t in SESS:
        E += [(s, nv), (t, nv+1)]; ss.append((nv, nv+1)); nv += 2
    return nv, E, ss
def minimal():
    """keep exactly one terminal at each terminal-hosting vertex (chosen by max matching vertex->session),
    move every other terminal to its own pendant leaf. Terminals stay pairwise distinct."""
    k = len(SESS); ends = defaultdict(list)
    for i, (q, s, t) in enumerate(SESS): ends[s].append((i, 0)); ends[t].append((i, 1))
    match_s = {}                                     # session -> vertex
    def aug(v, seen):
        for i, _ in ends[v]:
            if i in seen: continue
            seen.add(i)
            if i not in match_s or aug(match_s[i], seen): match_s[i] = v; return True
        return False
    for v in ends: aug(v, set())
    keep = {}                                        # vertex -> (session, side)
    for i, v in match_s.items():
        side = 0 if SESS[i][1] == v else 1; keep[v] = (i, side)
    for v in ends:                                   # vertices not matched still keep one terminal
        if v not in keep: keep[v] = ends[v][0]
    E = list(CORE); nv = 182; ss = []
    for i, (q, s, t) in enumerate(SESS):
        pair = []
        for side, x in ((0, s), (1, t)):
            if keep[x] == (i, side): pair.append(x)
            else: E.append((x, nv)); pair.append(nv); nv += 1
        ss.append(tuple(pair))
    allterm = [x for p in ss for x in p]; assert len(set(allterm)) == 2*k
    return nv, E, ss
if __name__ == "__main__":
    for f in (full, minimal):
        nv, E, ss = f()
        nl = len(E)-735
        from collections import Counter
        hop = Counter((s >= 182) + (t >= 182) for s, t in ss)
        print(f.__name__, "vertices", nv, "edges", len(E), "leaves", nl, "sessions by #leaves", sorted(hop.items()))
