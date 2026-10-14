"""Independent checker for Lemma B' certificates (sympy only; geometry rebuilt in check_C.py).

Lemma B': in a class block X = class(a) (live lines), P = U V^T (U, V in K^{X x 2}) plus at most one
rank-1 term u v^T with u, v supported on W_x = {r in X : r perp x}, the masked equations
   M[r,r] = 4,  M[r,s] = 0 (r perp s)
have no solution over any field of characteristic != 2 (cofactors lie in Z[1/2]; checked below).
Gauge: r0, r1 in X \\ W_x, r0 perp r1 -> U_r0, U_r1 independent (U_r0.V_r0 = 4, U_r1.V_r0 = 0), so
GL2 makes U_r0 = (1,0), U_r1 = (0,1).  Scaling of (u,v): u[w0] in {0,1}.
For each (case, class a, term x, b) the file certB/... gives pairs (f_i, g_i).  We check:
  * every f_i is one of the block equations (rebuilt here) or one of the gauge equations;
  * r0, r1, w0 satisfy the gauge hypotheses;
  * sum g_i f_i == 1 exactly over QQ.
"""
import os, sys, itertools
import sympy as sp
sys.argv = sys.argv[:1]
from check_C import L, ip, F0, LIVE, OPEN


def block(A, c, a, x):
    two = [q for q, k in zip(A, c) if k == 2]
    X = [r for r in LIVE if r not in two and ip(r, a) == 0 and any(ip(r, b) != 0 for b in two)]
    W = [r for r in X if ip(r, x) == 0]
    eqs = set()
    for r in X:
        for s in X:
            if r == s or ip(r, s) == 0:
                t = f'U{r}_0*V{s}_0 + U{r}_1*V{s}_1'
                if r in W and s in W:
                    t += f' + u{x}_{r}*v{x}_{s}'
                eqs.add(sp.expand(sp.sympify(t + (' - 4' if r == s else ''))))
    return X, W, eqs


def main():
    files = sorted(os.listdir('certB'))
    seen = set(); n = 0
    for ci, (A, c) in enumerate(OPEN):
        two = [q for q, k in zip(A, c) if k == 2]; one = [q for q, k in zip(A, c) if k == 1]
        for a in two:
            for x in one:
                X, W, eqs = block(A, c, a, x)
                for b in (0, 1):
                    fn = f'certB/{ci}_{a}_{x}_{b}.txt'
                    pairs = [ln.split(' ||| ') for ln in open(fn).read().strip().split('\n')]
                    total = 0; gauge_lines = set()
                    for f, g in pairs:
                        fe = sp.expand(sp.sympify(f)); ge = sp.sympify(g)
                        # cofactor denominators are powers of 2 -> identity holds over Z[1/2], any char != 2
                        for co in sp.expand(ge).as_coefficients_dict().values():
                            q = sp.Rational(co).q; assert q & (q - 1) == 0, (fn, co)
                        if fe not in eqs:
                            gauge_lines.add(f.strip())
                        total += ge * fe
                    # identify and validate gauge equations
                    r0s, r1s, w0s = set(), set(), set()
                    for gl in gauge_lines:
                        if gl.startswith('u'):
                            var, val = gl.split(' - ')
                            assert val == str(b) and var.split('_')[0] == f'u{x}'; w0s.add(int(var.split('_')[1]))
                        else:
                            var = gl.split(' ')[0]; idx = int(var.split('_')[0][1:])
                            if gl == f'U{idx}_0 - 1' or gl == f'U{idx}_1': r0s.add(idx)
                            elif gl == f'U{idx}_0' or gl == f'U{idx}_1 - 1': r1s.add(idx)
                            else: raise AssertionError(gl)
                    assert len(r0s) <= 1 and len(r1s) <= 1 and len(w0s) <= 1, fn
                    r0 = min(r0s) if r0s else None; r1 = min(r1s) if r1s else None; w0 = min(w0s) if w0s else None
                    for q in (r0, r1):
                        if q is not None:
                            assert q in X and q not in W   # P alone on its diagonal
                    if r0 is not None and r1 is not None:
                        assert ip(r0, r1) == 0 and r0 != r1
                    if w0 is not None:
                        assert w0 in W
                    assert sp.expand(total) == 1, fn
                    seen.add(fn); n += 1
    assert len(seen) == 414 == len(files)
    print('LEMMA B\' certificates verified:', n, '(all cofactor denominators are powers of 2)', '(23 cases x 3 classes x 3 terms x 2 scaling cases)')


main()
