"""W(E8)-orbits of pairwise-odd line sets (cliques) S, and the Clean(C(S)) test.
Clean(K): for every class O (h = 1..510): rank_Q R[O minus K] = 8, or rank 7 with kernel u ~ r_f, O = E7A1(f), f in K.
Clean(K) => the deficient lines of a kappa-side with support K are exactly e9 and the standard g_f (README Lemma 2);
Clean(C(S)) => Clean(B) for every B subset of C(S)."""
import numpy as np, json, sys
from e8 import *
from cliques_lib import NO, common
from catalog import SIMPLE
R = np.array(ROOT, float)
idx = {r: i for i, r in enumerate(ROOT)}
def refl(a):
    al = SIMPLE[a]; aa = sum(x * x for x in al)
    perm = []
    for r in ROOT:
        d = sum(x * y for x, y in zip(r, al)); s = tuple(x - 2 * d * y // aa for x, y in zip(r, al))
        perm.append(idx[s] if s in idx else idx[tuple(-x for x in s)])
    return perm
GEN = [refl(i) for i in range(8)]
for g in GEN: assert all(IP[g[a]][g[b]] == IP[a][b] or IP[g[a]][g[b]] == -IP[a][b] for a in range(0, NL, 7) for b in range(NL))
def apply(g, m): return sum(1 << g[a] for a in bits(m))
def cliques(k):
    out = []
    def rec(cl, cand):
        if len(cl) == k: out.append(sum(1 << c for c in cl)); return
        for c in bits(cand): rec(cl + [c], cand & NO[c] & ~((1 << (c + 1)) - 1))
    for a in range(NL): rec([a], NO[a] & ~((1 << (a + 1)) - 1))
    return out
def orbits(k):
    cl = cliques(k); pos = {m: i for i, m in enumerate(cl)}; par = list(range(len(cl)))
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for i, m in enumerate(cl):
        for g in GEN:
            j = pos[apply(g, m)]; a, b = find(i), find(j)
            if a != b: par[a] = b
    reps = {}
    for i in range(len(cl)): reps.setdefault(find(i), []).append(cl[i])
    return [(len(v), min(v)) for v in reps.values()]
CL = {h: sum(1 << a for a in range(NL) if dotq(LAB[a], h)) for h in range(1, 511)}
E7A1 = [ORTH[f] | (1 << f) for f in range(NL)]
def clean(K):
    bad = []
    for h, O in CL.items():
        L = list(bits(O & ~K)); M = R[L]
        rk = np.linalg.matrix_rank(M) if L else 0
        if rk == 8: continue
        if rk == 7:
            u = np.linalg.svd(M)[2][-1]
            f = [f for f in range(NL) if abs(abs(R[f] @ u) / np.linalg.norm(R[f]) - 1) < 1e-9]
            if f and O == E7A1[f[0]] and (K >> f[0]) & 1: continue
        bad.append((h, rk))
    return bad
if __name__ == '__main__':
    res = {}
    for k in range(1, 9):
        orb = orbits(k)
        for size, m in sorted(orb):
            S = list(bits(m)); C = common(S); Cm = sum(1 << c for c in C)
            b = clean(Cm)
            rkS = np.linalg.matrix_rank(R[S])
            print(f'k={k} orbit size {size} rep {S} |C(S)|={len(C)} rank(S)={rkS} Clean(C(S)) bad classes: {len(b)}', b[:4], flush=True)
            res.setdefault(k, []).append({'size': int(size), 'rep': [int(x) for x in S], 'C': [int(x) for x in C], 'rankS': int(rkS), 'bad': [(int(h), int(r)) for h, r in b]})
    json.dump(res, open('orbits.json', 'w'))
