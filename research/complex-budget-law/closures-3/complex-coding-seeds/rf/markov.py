"""Moment (Chebyshev-Markov) ceiling for char-0 Tanner gaps with half-rate flat-leverage (WH) local codes.
On a girth-2g layout, the per-edge spectral measure nu of Pi_P Pi_L Pi_P has moments j = 0..g-1 equal to the
free-projection law nu_free = 1/2 delta_0 + 1/2 arcsine[0,1]  (closed walks of length < 2g are tree-like).
S0/E = nu({1}) <= M_g := min { int p dnu_free : deg p <= g-1, p >= 0 on [0,1], p(1) = 1 }.
Threshold for a > 0 with b = girth - 1 = 2g - 1:  S/E > 1/(2g).  Prints M_g, 1/(2g), and a_max = 1 - ln((1-M)/M)/ln(2g-1)."""
import numpy as np
from math import comb, log
from scipy.optimize import linprog

xs = np.linspace(0, 1, 4001)
ok = True
for g in range(2, 13):
    d = g - 1
    mom = np.array([1.0] + [0.5 * comb(2 * k, k) / 4 ** k for k in range(1, d + 1)])  # int x^k dnu_free
    # p(x) = sum c_k x^k  (monomial basis is fine for d <= 11 on [0,1] with exact LP)
    V = np.vstack([xs ** k for k in range(d + 1)]).T
    res = linprog(mom, A_ub=-V, b_ub=np.zeros(len(xs)), A_eq=np.ones((1, d + 1)), b_eq=[1.0],
                  bounds=[(None, None)] * (d + 1), method='highs')
    M = res.fun; thr = 1 / (2 * g)
    a = 1 - log((1 - M) / M) / log(2 * g - 1)
    ok &= (abs(M - thr) < 1e-3) if g % 2 == 0 else (M > thr and (g != 3 or abs(M - 2/11) < 1e-3))
    print(f'girth {2*g:2d}  b={2*g-1:2d}  M_g={M:.4f}  thr=1/(2g)={thr:.4f}  ratio={M/thr:.3f}  a_max={a:+.3f}')
print('PASS_EXPRF_MARKOV' if ok else 'FAIL_EXPRF_MARKOV')   # M_g = 1/(2g) for even g (a <= 0); girth 6: 2/11
