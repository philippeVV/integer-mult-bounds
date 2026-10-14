"""Exact checks of the 2-adic Smith lemma (smith.py)."""
import sys, json, random, os
from smith import *
from common import *
HERE = os.path.dirname(os.path.abspath(__file__))
def Bt_block(basis, R, A):
    return [[dot(x, y) for b in A for y in basis[b]] for a in R for x in basis[a]]
def lemma_sum(basis, terms):
    tot = 0; rows = []
    for d, t in terms:
        e = min(v2(x) for x in d.values()) + min(v2(y) for y in t.values())
        R, A = sorted(d), sorted(t)
        M = Bt_block(basis, R, A)
        rho = r_j(M, 0); dd = max(0, -e)
        h = r_j(M, dd) if e > -10 else None
        rows.append((len(R), len(A), rho, e, h, g(dd, rho)))
        assert h <= g(dd, rho)
        tot += h
    return tot, rows
# (0) sanity: r_j against brute-force definition on tiny random integer matrices: r_j(X+Y) <= r_j(X)+r_j(Y), r_0 = F2 rank
rng = random.Random(1)
for _ in range(300):
    a, b = rng.randint(1, 6), rng.randint(1, 6)
    X = [[rng.randint(-9, 9) for _ in range(b)] for _ in range(a)]; Y = [[rng.randint(-9, 9) for _ in range(b)] for _ in range(a)]
    S = [[X[i][j] + Y[i][j] for j in range(b)] for i in range(a)]
    for j in range(4): assert r_j(S, j) <= r_j(X, j) + r_j(Y, j)
    assert r_j(X, 0) == rank2([sum(((X[i][k] & 1) << k) for k in range(b)) for i in range(a)])
print('PASS r_j subadditivity / r_0 = F2 rank on 300 random pairs')
# (1) g_d bound on random rectangles of random plane assignments
RC = [rchoices(q) for q in P]; worst = {}
for _ in range(400):
    r = [rng.choice(RC[i]) for i in range(n)]; basis = [(P[i], r[i]) for i in range(n)]
    R = rng.sample(range(n), rng.randint(1, 30)); A = rng.sample(range(n), rng.randint(1, 30))
    M = Bt_block(basis, R, A); rho = r_j(M, 0)
    for d in range(3):
        h = r_j(M, d); assert h <= g(d, rho), (h, d, rho)
        worst[(rho, d)] = max(worst.get((rho, d), 0), h)
print('PASS r_d(Bt[R,A]) <= g_d(rho) on 400 random rectangles; max seen (rho,d)->h:', dict(sorted(worst.items())))
# (2) naive decoder: equality 224
r = [RC[i][0] for i in range(n)]; basis = [(P[i], r[i]) for i in range(n)]
tot, _ = lemma_sum(basis, [({i: 1}, {i: 1}) for i in range(n)]); print('naive decoder: sum r =', tot, '(lemma: >= 224)')
# (3) exact plane decoders from expJ (rect*/rectF*), verified exactly first
for f in sys.argv[1:]:
    dd = json.load(open(f)); r = dd['r']; basis = [(P[i], r[i]) for i in range(n)]
    terms = [({int(i): Fr(x) for i, x in a.items()}, {int(j): Fr(y) for j, y in b.items()}) for a, b in dd['terms']]
    bad, C, _ = verify(r, terms); assert bad == 0
    tot, rows = lemma_sum(basis, terms)
    print(os.path.basename(f), 'exact OK, C =', C, ' Smith sum =', tot, '>= 224:', tot >= 224, ' depths', sorted(set(-x[3] for x in rows)), flush=True)
