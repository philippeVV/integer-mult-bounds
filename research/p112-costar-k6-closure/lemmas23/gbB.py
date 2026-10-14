"""Lemma B': class block with P (rank 2) + at most one rank-1 term Z_x (supported on W_x) is infeasible.
Groebner basis of the block system (only term x kept) must be {1}."""
import sys, time
from sage.all__sagemath_singular import PolynomialRing, QQ, GF
from cases import *


def system(A, c, a, keep):
    two, one = split(A, c); E6, cl = classes(two)
    X = [r for r in cl[a] if r not in F0]
    names = []; var = {}
    for r in X:
        for k in range(2):
            var['U', r, k] = f'U{r}_{k}'; var['V', r, k] = f'V{r}_{k}'; names += [var['U', r, k], var['V', r, k]]
    for x in keep:
        for r in X:
            if GM[r][x] == 0:
                var['u', x, r] = f'u{x}_{r}'; var['v', x, r] = f'v{x}_{r}'; names += [var['u', x, r], var['v', x, r]]
    eqs = []
    for r in X:
        for s in X:
            if r == s or GM[r][s] == 0:
                t = [f'{var["U", r, k]}*{var["V", s, k]}' for k in range(2)]
                t += [f'{var["u", x, r]}*{var["v", x, s]}' for x in keep if ('u', x, r) in var and ('v', x, s) in var]
                eqs.append(' + '.join(t) + (' - 4' if r == s else ''))
    return X, names, eqs


if __name__ == '__main__':
    ci = int(sys.argv[1]); p = int(sys.argv[2])
    A, c = OPEN[ci]; two, one = split(A, c)
    for a in two:
        for x in one:
            X, names, eqs = system(A, c, a, [x])
            R = PolynomialRing(QQ if p == 0 else GF(p), names, order='degrevlex')
            t = time.time(); gb = R.ideal([R(e) for e in eqs]).groebner_basis()
            print(ci, a, x, len(X), len(names), len(eqs), 'p', p,
                  'ONE' if list(gb) == [1] else 'NOT_ONE len %d' % len(gb), '%.1fs' % (time.time() - t), flush=True)
