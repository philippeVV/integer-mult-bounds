"""Orbits of covering normal sets under the frame stabilizer H (order 172,032)."""
import pickle, sys
from collections import defaultdict
from model import F0

GENS = pickle.load(open("gens_small.pkl", "rb"))

def orbits_of_sets(sets):
    idx = {s: i for i, s in enumerate(sets)}
    par = list(range(len(sets)))
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for i, s in enumerate(sets):
        for g in GENS:
            t = tuple(sorted(g[a] for a in s))
            j = idx[t]   # KeyError would mean the family is not H-invariant
            ri, rj = find(i), find(j)
            if ri != rj:
                par[ri] = rj
    orb = defaultdict(list)
    for i in range(len(sets)):
        orb[find(i)].append(i)
    out = [(sets[min(v)], len(v), [sets[i] for i in v]) for v in orb.values()]
    out.sort(key=lambda x: x[0])
    return out

if __name__ == "__main__":
    from cover import covering_sets
    c3 = covering_sets(3)
    o3 = orbits_of_sets(c3)
    print("covering 3-set orbits:", len(o3), [(r, n, sum(1 for a in r if a in F0)) for r, n, _ in o3])
    c4 = pickle.load(open("cover4.pkl", "rb"))
    o4 = orbits_of_sets(c4)
    print("covering 4-set orbits:", len(o4))
    print("sizes:", sorted(n for _, n, _ in o4))
    pickle.dump([(r, n) for r, n, _ in o3], open("orb3.pkl", "wb"))
    pickle.dump([(r, n) for r, n, _ in o4], open("orb4.pkl", "wb"))
