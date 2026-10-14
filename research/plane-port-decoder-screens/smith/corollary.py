"""Universal (every plane assignment) class bounds from the Smith lemma:  224 <= sum_k h_k,  h_k <= cap(rho_k, depth_k),
cap = min(g_d(rho), exact H(rho,d) for rho in {8,9}),  depth 0 (2-integral) -> cap = rho.  Hence C >= 224 / max cap/rho."""
from fractions import Fraction as Fr
from math import ceil
from smith import g
import hcap   # prints exact H for rho=9,8 (d=0..3)
H = {9: hcap.out[(9, 'V,V')], 8: [max(hcap.out[k][d] for k in hcap.out if k[0] == 8) for d in range(4)]}
def cap(rho, d):
    c = g(d, rho) if d <= 3 else 2 ** rho - 1
    if rho in H and d <= 3: c = min(c, H[rho][d])
    return c
def bound(rhos, depths):
    worst = max(Fr(cap(r, d), r) for r in rhos for d in depths)
    return worst, ceil(Fr(224) / worst - Fr(1, 10**9))
print('cap table (rho: d=0..3):', {r: [cap(r, d) for d in range(4)] for r in range(1, 10)})
w, b = bound(range(1, 10), [0, 1]); print('all terms half-integral (depth<=1):             max cap/rho =', w, '=> C >=', b)
w, b = bound(range(1, 10), [0, 1, 2]); print('depth<=2:                                     max cap/rho =', float(w), '=> C >=', b)
for k in (2, 3, 4, 5):
    w = max(Fr(2 ** r - 1, r) for r in range(1, k + 1)); print(f'every non-2-integral term has rho<={k} (any depth): max cap/rho = {w} => C >= {ceil(Fr(224)/w)}')
w, b = bound(range(1, 6), [0, 1]); print('depth<=1 and rho<=5: C >=', b, '; depth<=1 and rho<=3: C >=', bound(range(1, 4), [0, 1])[1])
