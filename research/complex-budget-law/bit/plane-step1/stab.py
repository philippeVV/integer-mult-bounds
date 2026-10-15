"""Stabilizer of the puncture P in O(q) on E8/2E8 (= W(E8)/{+-1} acting on lines). Output: permutations of 120 lines."""
import core, json
P = list(range(111, 120)); V = core.VEC
print('rank P', core.rank2([V[t] for t in P]))
# choose independent subset basis of P
basis = []
for t in P:
    if core.rank2([V[b] for b in basis] + [V[t]]) > len(basis): basis.append(t)
print('basis', basis)
def solve_coords(v, bas):  # express v in span of bas vectors
    import itertools
    for m in range(1 << len(bas)):
        x = 0
        for i, b in enumerate(bas):
            if m >> i & 1: x ^= V[b]
        if x == v: return m
    return None
# need full basis of F2^8 to define linear map; if P has rank 8 fine
assert len(basis) == 8
coord = [solve_coords(V[t], basis) for t in range(120)]
perms = []
def extend(img):
    if len(img) == len(P):
        mp = dict(zip(P, img))
        # linear map from basis images
        def f(m):
            x = 0
            for i, b in enumerate(basis):
                if m >> i & 1: x ^= V[mp[b]]
            return x
        perm = []
        for t in range(120):
            u = f(coord[t])
            if u not in core.VIDX: return
            perm.append(core.VIDX[u])
        if sorted(perm) != list(range(120)): return
        if sorted(perm[t] for t in P) != P: return
        perms.append(perm); return
    i = len(img); a = P[i]
    for c in P:
        if c in img: continue
        if all(core.G[a][P[j]] % 2 == core.G[c][img[j]] % 2 for j in range(i)):
            extend(img + [c])
extend([])
# check isometry: B preserved on all lines
for perm in perms:
    assert all(core.G[s][t] % 2 == core.G[perm[s]][perm[t]] % 2 for s in range(120) for t in range(0, 120, 7))
print('stabilizer order', len(perms))
json.dump(perms, open('stabP.json', 'w'))
