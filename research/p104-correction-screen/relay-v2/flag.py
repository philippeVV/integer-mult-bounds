"""Global-flag decision-tree designs: evaluate exact-ish junk J for an order of 8 functionals."""
import json, random, sys, collections
import numpy as np
D = json.load(open('geo.json')); n = D['n']; q = D['q']; nt = D['Hnormal']; G = np.array(D['G'])
C = np.eye(n) - G / 2
dot = lambda a, b: bin(a & b).count('1') & 1
US = [u for u in range(1, 512) if all(dot(u, x) for x in q)][0]
supp = [[s for s in range(n) if not dot(q[s], nt[t])] for t in range(n)]
def rank2(vs):
    b = []
    for v in vs:
        for x in b: v = min(v, v ^ x)
        if v: b.append(v)
    return len(b)
def evaluate(zs, detail=False):
    code = [[dot(z, q[s]) for z in zs] for s in range(n)]
    users = collections.defaultdict(set); kids = {}
    def walk(t, S, k):
        # S: tuple of sources (len>=2); find first level >=k that splits S
        while k < len(zs):
            a = tuple(s for s in S if code[s][k] == 0); b = tuple(s for s in S if code[s][k] == 1)
            if a and b: break
            k += 1
        assert k < len(zs)
        users[S].add(t); kids[S] = (a, b)
        for c in (a, b):
            if len(c) >= 2: walk(t, c, k + 1)
    tops = {}
    for t in range(n):
        S = tuple(supp[t]); k = 0
        while True:
            a = tuple(s for s in S if code[s][k] == 0); b = tuple(s for s in S if code[s][k] == 1)
            if a and b: break
            k += 1
        tops[t] = (a, b)
        for c in (a, b):
            if len(c) >= 2: walk(t, c, k + 1)
    J = n  # targets: A_t - 1 = 1 each
    per = collections.Counter()
    for S, us in users.items():
        M = C[np.ix_(sorted(us), S)]
        r = np.linalg.matrix_rank(M)
        a, b = kids[S]; ia = [S.index(s) for s in a]; ib = [S.index(s) for s in b]
        ra = np.linalg.matrix_rank(M[:, ia]); rb = np.linalg.matrix_rank(M[:, ib])
        J += ra + rb - r; per[len(S)] += ra + rb - r
    if detail: return J, users, kids, tops, per
    return J
def random_basis(rng):
    while True:
        zs = rng.sample(range(1, 512), 8)
        if rank2(zs + [US]) == 9: return zs
if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1])); iters = int(sys.argv[2])
    best = None
    zs = random_basis(rng); cur = evaluate(zs)
    for it in range(iters):
        z2 = list(zs)
        mv = rng.random()
        if mv < 0.4:
            i, j = rng.sample(range(8), 2); z2[i], z2[j] = z2[j], z2[i]
        else:
            i = rng.randrange(8); z2[i] = rng.randrange(1, 512)
            if rank2(z2 + [US]) < 9: continue
        v = evaluate(z2)
        if v <= cur or rng.random() < 0.02:
            zs, cur = z2, v
        if best is None or cur < best[0]:
            best = (cur, list(zs)); print(it, cur, flush=True)
    print('BEST', best)
