"""Abstract triangle system T (independent of E8) and a Nullstellensatz certificate sum g_i f_i = 1 over QQ.
Cliques Q1, Q3, Q5, each 4 rows split into two parts; every row r has u_r, v_r in K^2 (the P-coordinates).
 Q1, Q3 (pure-like):  sum_{r in Q} v_r u_r^T = 4 I_2                      (4 eqs each)
 Q5 (deficient):      det(4 I_2 - sum_{r in Q5} v_r u_r^T) = 0             (1 eq)
 cross, matched parts (A,C),(A',C'),(A,B),(A',B'),(C,B),(C',B'), both directions: u_r . v_s = 0.
Output: certificate pickle with polynomials as {exponent-tuple: Fraction} dicts."""
import pickle, time, sys
from fractions import Fraction
from sage.all__sagemath_singular import PolynomialRing, QQ

ROWS = {'1': ['a1', 'a2', 'a3', 'a4'], '3': ['c1', 'c2', 'c3', 'c4'], '5': ['b1', 'b2', 'b3', 'b4']}
PARTS = {'1': (['a1', 'a2'], ['a3', 'a4']), '3': (['c1', 'c2'], ['c3', 'c4']), '5': (['b1', 'b2'], ['b3', 'b4'])}
names = [f'{t}{r}_{k}' for q in '135' for r in ROWS[q] for t in 'uv' for k in range(2)]

def equations(R):
    u = lambda r, k: R(f'u{r}_{k}'); v = lambda r, k: R(f'v{r}_{k}')
    E = []
    for q in '13':
        for g in range(2):
            for h in range(2):
                E.append(sum(v(r, g) * u(r, h) for r in ROWS[q]) - (4 if g == h else 0))
    N = [[sum(v(r, g) * u(r, h) for r in ROWS['5']) for h in range(2)] for g in range(2)]
    E.append((4 - N[0][0]) * (4 - N[1][1]) - N[0][1] * N[1][0])
    for x, y in (('1', '3'), ('1', '5'), ('3', '5')):
        for i in range(2):
            for r in PARTS[x][i]:
                for s in PARTS[y][i]:
                    E.append(sum(u(r, k) * v(s, k) for k in range(2)))
                    E.append(sum(u(s, k) * v(r, k) for k in range(2)))
    return E

if __name__ == '__main__':
    R = PolynomialRing(QQ, names, order='degrevlex')
    E = equations(R)
    t = time.time()
    I = R.ideal(E)
    G = I.groebner_basis(); print('GB', list(G) == [1], '%.1fs' % (time.time() - t), flush=True)
    t = time.time()
    co = R(1).lift(I)
    print('lift %.1fs' % (time.time() - t), flush=True)
    assert sum(c * f for c, f in zip(co, E)) == 1
    dens = set()
    out = []
    for c in co:
        d = {}
        for e, cf in c.dict().items():
            q = Fraction(int(cf.numerator()), int(cf.denominator())); d[tuple(e)] = q; dens.add(q.denominator)
        out.append(d)
    pw2 = all(dd & (dd - 1) == 0 for dd in dens)
    print('terms', sum(len(d) for d in out), 'max deg', max((sum(e) for d in out for e in d), default=0),
          'denominators powers of 2:', pw2, sorted(dens)[-5:])
    pickle.dump({'names': names, 'cofactors': out}, open('triangle_cert.pkl', 'wb'))
