"""Lemma B' (gauge-fixed): a class block X with P = U V^T (rank <= 2) plus at most ONE rank-1 term
Z_x = u v^T (u, v supported on W_x = {r in X : r perp x}) cannot give 4 on the diagonal and 0 on
orthogonal pairs.

Gauge (no loss of generality, no case split): pick r0, r1 in N_x = X \\ W_x with r0 perp r1.
Their diagonals and the pair (r0,r1),(r1,r0) receive only P, so U_r0.V_r0 = 4, U_r1.V_r1 = 4,
U_r0.V_r1 = U_r1.V_r0 = 0.  Hence U_r0, U_r1 are linearly independent and by GL2 we may take
U_r0 = (1,0), U_r1 = (0,1) (then V_r0 = (4,0), V_r1 = (0,4) follow from the equations).
Rank-1 scaling u -> t u, v -> v / t: split u[w0] = 1 or u[w0] = 0 (exhaustive).
Each case: Groebner basis == {1}  ==> no solution over any field containing QQ (p = 0)."""
import sys, time
from sage.all__sagemath_singular import PolynomialRing, QQ, GF
from cases import *
from gbB import system


def pick(A, c, a, x):
    two, one = split(A, c); E6, cl = classes(two)
    X = [r for r in cl[a] if r not in F0]
    Nx = [r for r in X if GM[r][x] != 0]
    for r0 in Nx:
        for r1 in Nx:
            if r0 < r1 and GM[r0][r1] == 0:
                return X, r0, r1, [r for r in X if GM[r][x] == 0][0]
    raise SystemExit('no orthogonal pair in N_x')


def run(ci, a, x, p):
    A, c = OPEN[ci]
    X, r0, r1, w0 = pick(A, c, a, x)
    X, names, eqs = system(A, c, a, [x])
    R = PolynomialRing(QQ if p == 0 else GF(p), names, order='degrevlex')
    base = [R(e) for e in eqs]
    g = [R(f'U{r0}_0') - 1, R(f'U{r0}_1'), R(f'U{r1}_0'), R(f'U{r1}_1') - 1]
    res = []
    for b in (1, 0):
        t = time.time()
        gb = R.ideal(base + g + [R(f'u{x}_{w0}') - b]).groebner_basis()
        res.append((b, list(gb) == [1], len(gb), time.time() - t))
    return (r0, r1, w0), res


if __name__ == '__main__':
    p = int(sys.argv[1]); cis = [int(v) for v in sys.argv[2:]] or range(len(OPEN))
    for ci in cis:
        A, c = OPEN[ci]; two, one = split(A, c)
        for a in two:
            for x in one:
                gauge, res = run(ci, a, x, p)
                ok = all(r[1] for r in res)
                print(ci, 'class', a, 'term', x, 'gauge', gauge, 'p', p, 'ONE' if ok else 'FAIL',
                      ['u=%d:%s/%.1fs' % (b, o, tt) for b, o, l, tt in res], flush=True)
