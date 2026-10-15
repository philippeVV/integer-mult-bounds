"""Non-uniform coverage law: tau = sum_a dim E_a unpaid port-coordinates (112 ports, 8/9 families per port per call).
Exact break-even centre charge with ceiling.py (word unchanged, favourable), and the dyadic gain a decoder needs
above the 2-integral floor C >= tau (Theorem A generalised: sum rho >= rank blockdiag Gram(E_a) = tau)."""
import json
from fractions import Fraction as Q
from common import *
base = cap(WORD)
def be(u):
    lo, hi = 0, 600
    while hi - lo > 1:
        mid = (lo + hi) // 2
        try: ok = cap(WORD, centre=mid, unpaid=u) > base
        except Exception: ok = False
        lo, hi = (mid, hi) if ok else (lo, mid)
    return lo
rows = []
for tau in (112, 120, 128, 144, 168, 224, 280, 336, 448):
    u = Q(tau, 112)
    b = be(u)
    rows.append(dict(tau=tau, unpaid_avg=str(u), break_even_C=b, first_order=float(64 + Q(8, 9) * (tau - 112)),
                     gain_needed=tau - b, gain_needed_per_coord=round((tau - b) / tau, 4)))
    print(rows[-1])
# kappa caps at selected favourable charges (complex side only; bit side similar x, binding)
caps = {}
for u, Cs in ((2, (25, 50, 96, 112, 127, 160)), (3, (46, 126, 200, 227, 260))):
    for Cc in Cs:
        a = cap(WORD, centre=Cc, unpaid=u)
        caps['u%d_C%d' % (u, Cc)] = dict(a=float(a), ratio=float(a / base), kappa=float(a / (1 + a)))
        print('u', u, 'C', Cc, 'a_cap %.6f x%.2f kappa_cap %.6f' % (float(a), float(a / base), float(a / (1 + a))))
json.dump(dict(rows=rows, caps=caps), open('d_price.json', 'w'), indent=1)
