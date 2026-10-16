#!/usr/bin/env python3
"""expNE: small-unit pair kappa when the s = 3 word realises eta3 (transfer.py) instead of 'cal' eta.
Unit found in complex orientation at window h (complex h, bit twin h-1, same x), as expNB table (a).
's3(rho)': side saving = eta3(h_side, x, rho) * X(3, h_side, x), with E8's calibration offset kept
           (bit/cx realised s=3 E8 ratios 0.1319/0.1350 vs transfer model 0.1504/0.141 at rho_today -> scale by those).
'cal9'   : expNB (word stays at s = 9 with hosting).  'best' = max(s3, cal9) per side.  'cal' = expNB 'cal'."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # kappa_pair.py (expNB), transfer.py (expNE)
import kappa_pair as KP
from transfer import eta3

RHO_T = {'bit': 3512 / 1554, 'cx': 4142 / 1792}
OFF = {'bit': 0.1319 / eta3(8, 63 / 111, RHO_T['bit']), 'cx': 0.1350 / eta3(9, 64 / 112, RHO_T['cx'])}


def side(sd, h, x, model, rho_scale=1.0):
    if model == 'cal9' or model == 'cal':
        return KP.side_saving(sd, h, x, model)
    s3 = OFF[sd] * eta3(h, x, RHO_T[sd] * rho_scale) * KP.X(3, h, x)
    if model == 's3':
        return s3
    return max(s3, KP.side_saving(sd, h, x, 'cal9'))


def pair(h, x, model, rs=1.0):
    a, b = side('bit', h - 1, x, model, rs), side('cx', h, x, model, rs)
    A = min(a, b)
    return A / (1 + A)


print('offsets (exact E8 s=3 / transfer model):', {k: round(v, 3) for k, v in OFF.items()})
print('E8 itself: s3 pair kappa %.5f (x%.2f); cal %.5f' % (pair(9, 64 / 112, 's3'), pair(9, 64 / 112, 's3') / 0.0017,
                                                          pair(9, 64 / 112, 'cal')))
xs = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55]
for title, model, rs in [('cal (expNB, eta 0.29/0.33 at best s)', 'cal', 1), ('cal9 (word kept at s=9, hosted)', 'cal9', 1),
                         ('s3, today helper burden rho = 2.26/2.31', 's3', 1),
                         ('best(s3 today, cal9)', 'best', 1),
                         ('s3, helper mass halved (rho ~1.1, the bit design-free floor 212/439)', 's3', 0.483),
                         ('best(s3 halved, cal9)', 'best', 0.483)]:
    print('\n[%s]  multiplier over 0.0017' % title)
    print('h_c  ' + ''.join('x=%-6.2f' % x for x in xs))
    for h in (5, 6, 7):
        print('%d    ' % h + ''.join('%-8.2f' % (pair(h, x, model, rs) / 0.0017) for x in xs))

# s = 5 with optimistic hosting (s5_bound.py, slot packing): eta5 = 0.2635 bit / 0.2613 cx of X(5); 'perfect' 0.2894 / 0.2870
ETA5 = {'slot': {'bit': 0.2635, 'cx': 0.2613}, 'perfect': {'bit': 0.2894, 'cx': 0.2870}}


def best_all(h, x, pack):
    def sd(s_, hh):
        return max(side(s_, hh, x, 'best'), ETA5[pack][s_] * KP.X(5, hh, x))
    A = min(sd('bit', h - 1), sd('cx', h))
    return A / (1 + A)


for pack in ('slot', 'perfect'):
    print('\n[best over s in {3 (realised-transfer), 5 (eta5 %s hosting), 9 (cal9)}]  multiplier over 0.0017' % pack)
    print('h_c  ' + ''.join('x=%-6.2f' % x for x in xs))
    for h in (5, 6, 7):
        print('%d    ' % h + ''.join('%-8.2f' % (best_all(h, x, pack) / 0.0017) for x in xs))
