"""Produce Nullstellensatz certificates for Lemma B': cofactors g_i (over QQ) with sum g_i f_i = 1,
where f_i are the gauge-fixed class-block equations (see gbB2.py).  Written to certB/<ci>_<a>_<x>_<b>.txt
as lines 'equation_string ||| cofactor_string'."""
import sys, os, time
from sage.all__sagemath_singular import PolynomialRing, QQ, Matrix
from cases import *
from gbB import system
from gbB2 import pick

os.makedirs('certB', exist_ok=True)
cis = [int(v) for v in sys.argv[1:]] or range(len(OPEN))
for ci in cis:
    A, c = OPEN[ci]; two, one = split(A, c)
    for a in two:
        for x in one:
            X, r0, r1, w0 = pick(A, c, a, x)
            X, names, eqs = system(A, c, a, [x])
            R = PolynomialRing(QQ, names, order='degrevlex')
            gauge = [f'U{r0}_0 - 1', f'U{r0}_1', f'U{r1}_0', f'U{r1}_1 - 1']
            for b in (1, 0):
                fn = f'certB/{ci}_{a}_{x}_{b}.txt'
                if os.path.exists(fn):
                    continue
                gens = eqs + gauge + [f'u{x}_{w0} - {b}']
                I = R.ideal([R(e) for e in gens])
                t = time.time()
                co = R.ideal([R(1)]).gens()[0].lift(I)
                chk = sum(g * R(f) for g, f in zip(co, gens))
                assert chk == 1
                with open(fn + '.tmp', 'w') as fh:
                    for f, g in zip(gens, co):
                        if g != 0:
                            fh.write(f + ' ||| ' + str(g) + '\n')
                os.rename(fn + '.tmp', fn)
                print(ci, a, x, b, 'lift ok', sum(1 for g in co if g != 0), 'nonzero cofactors', 'max terms', max(len(g.monomials()) for g in co), '%.1fs' % (time.time() - t), flush=True)
