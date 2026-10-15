"""Cell-tree designs: a binary decision tree of functionals over source cells; target t's pools are supp(t) ∩ cell.
Evaluates junk J = 104 + sum_pools (rank(D|a) + rank(D|b) - rank(D))."""
import json, random, sys, collections, pickle
import numpy as np
D = json.load(open('geo.json')); n = D['n']; q = D['q']; nt = D['Hnormal']; G = np.array(D['G'])
C = np.eye(n) - G / 2
dot = lambda a, b: bin(a & b).count('1') & 1
supp = [tuple(s for s in range(n) if not dot(q[s], nt[t])) for t in range(n)]
ALL = tuple(range(n))
class Tree:
    def __init__(self, glist, zmap=None):
        self.glist = glist; self.zmap = dict(zmap or {})
    def z(self, cell):
        if cell in self.zmap: return self.zmap[cell]
        for z in self.glist:
            a = [s for s in cell if dot(z, q[s]) == 0]
            if a and len(a) < len(cell): return z
        for z in range(1, 512):
            a = [s for s in cell if dot(z, q[s]) == 0]
            if a and len(a) < len(cell): return z
    def split(self, cell):
        z = self.z(cell)
        a = tuple(s for s in cell if dot(z, q[s]) == 0); b = tuple(s for s in cell if dot(z, q[s]) == 1)
        return a, b
rcache = {}
def rk(M):
    return int(np.linalg.matrix_rank(M)) if M.size else 0
def evaluate(tr, detail=False):
    users = collections.defaultdict(set); kids = {}; cells_used = collections.Counter(); scache = {}
    def sp(cell):
        if cell not in scache: scache[cell] = tr.split(cell)
        return scache[cell]
    def walk(t, S, cell, top=False):
        while True:
            ca, cb = sp(cell)
            sa = set(ca); a = tuple(s for s in S if s in sa); b = tuple(s for s in S if s not in sa)
            if a and b: break
            cell = ca if a else cb
        cells_used[cell] += 1
        if top: kids[('H', t)] = (a, b)
        else: users[S].add(t); kids[S] = (a, b)
        for c, cc in ((a, ca), (b, cb)):
            if len(c) >= 2: walk(t, c, cc)
    for t in range(n): walk(t, supp[t], ALL, top=True)
    J = n; per = collections.Counter()
    for S, us in users.items():
        M = C[np.ix_(sorted(us), S)]; a, b = kids[S]
        ia = [S.index(s) for s in a]; ib = [S.index(s) for s in b]
        j = rk(M[:, ia]) + rk(M[:, ib]) - rk(M); J += j; per[len(S)] += j
    if detail: return J, users, kids, per, cells_used
    return J, cells_used
def search(tr, rng, iters, log):
    cur, cu = evaluate(tr); best = (cur, dict(tr.zmap))
    for it in range(iters):
        cells = list(cu)
        cell = rng.choice(cells)
        old = tr.zmap.get(cell, None)
        z = rng.randrange(1, 512)
        a = [s for s in cell if dot(z, q[s]) == 0]
        if not a or len(a) == len(cell): continue
        tr.zmap[cell] = z
        v, cu2 = evaluate(tr)
        if v <= cur or rng.random() < 0.01:
            cur, cu = v, cu2
            if v < best[0]: best = (v, dict(tr.zmap)); print(it, v, flush=True, file=log)
        else:
            if old is None: del tr.zmap[cell]
            else: tr.zmap[cell] = old
    return best
if __name__ == '__main__':
    seed = int(sys.argv[1]); iters = int(sys.argv[2]); glist = eval(sys.argv[3]); out = sys.argv[4]
    rng = random.Random(seed); tr = Tree(glist)
    if len(sys.argv) > 5:
        tr.zmap = pickle.load(open(sys.argv[5], 'rb'))[1]
    best = search(tr, rng, iters, sys.stdout)
    pickle.dump((best[0], best[1], glist), open(out, 'wb')); print('BEST', best[0])
