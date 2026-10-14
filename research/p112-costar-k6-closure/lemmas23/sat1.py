"""Support-level SAT relaxation on the union of the three A2-class blocks.
Booleans: U[x,r] = (u_x[r] != 0), V[x,r] = (v_x[r] != 0) for r in W_x (live class lines).
(1) cross-class orthogonal pair (r,s): Q[r,s] = sum_{x in O(r)&O(s)} u_x[r] v_x[s] = 0
    => the number of x with U[x,r]&V[x,s] is never exactly 1.
(2) every orthogonal triangle T inside one class: P_T = 4I - Q_T has rank <= 2, so det(4I - Q_T) = 0,
    so Q_T is not nilpotent-support: some cycle (loop, 2-cycle or 3-cycle) of possibly-nonzero Q entries.
    Q[r,s] possibly nonzero only if some x has U[x,r]&V[x,s]."""
import sys, itertools
from pycryptosat import Solver
from cases import *
def build(ci, extra=()):
    A, c = OPEN[ci]; two, one = split(A, c); E6, cl = classes(two)
    cls = {r: a for a in two for r in cl[a] if r not in F0}
    U = list(cls)
    O = {r: [x for x in one if GM[r][x] == 0] for r in U}
    nv = [0]; var = {}
    def new():
        nv[0] += 1; return nv[0]
    for r in U:
        for x in O[r]:
            var['U', x, r] = new(); var['V', x, r] = new()
    cl_ = []
    def prod(x, r, s):  # aux p <-> U[x,r] & V[x,s]
        k = ('p', x, r, s)
        if k not in var:
            p = new(); var[k] = p; a_, b_ = var['U', x, r], var['V', x, s]
            cl_.extend([[-p, a_], [-p, b_], [p, -a_, -b_]])
        return var[k]
    def poss(r, s):  # aux q -> OR_x p(x,r,s)   (q true only if some product nonzero)
        k = ('q', r, s)
        if k not in var:
            q = new(); var[k] = q
            cl_.append([-q] + [prod(x, r, s) for x in O[r] if x in O[s]])
        return var[k]
    ncross = 0
    for r in U:
        for s in U:
            if cls[r] != cls[s] and GM[r][s] == 0:
                xs = [x for x in O[r] if x in O[s]]
                ps = [prod(x, r, s) for x in xs]
                for i in range(len(ps)):   # not exactly one true
                    cl_.append([-ps[i]] + [ps[j] for j in range(len(ps)) if j != i])
                ncross += 1
    ntri = 0
    for a in two:
        L = [r for r in U if cls[r] == a]
        for T in itertools.combinations(L, 3):
            if all(GM[p][q] == 0 for p, q in itertools.combinations(T, 2)):
                ntri += 1
                cyc = []
                for r in T: cyc.append([poss(r, r)])
                for r, s in itertools.combinations(T, 2): cyc.append([poss(r, s), poss(s, r)])
                for r, s, t in itertools.permutations(T):
                    if r == min(T): cyc.append([poss(r, s), poss(s, t), poss(t, r)])
                # OR of conjunctions
                ors = []
                for conj in cyc:
                    z = new(); ors.append(z)
                    for l in conj: cl_.append([-z, l])
                cl_.append(ors)
    return cl_, nv[0], ncross, ntri, var
if __name__ == '__main__':
    for ci in range(len(OPEN)):
        cls_, n, nc, nt, var = build(ci)
        s = Solver(); s.add_clauses(cls_)
        sat, sol = s.solve()
        print(ci, OPEN[ci], 'vars', n, 'cross', nc, 'triangles', nt, 'SAT' if sat else 'UNSAT', flush=True)
