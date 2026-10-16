#!/usr/bin/env python3
"""expL2 rho_min.py: what the coupling lemma forces on rho = helper mass / data mass (per stage, data = 2*V*a, a = h-1),
and what the resulting small-unit pair multiplier is (expNE s=3 transfer, expNB pair; 'best' = max(s3, cal9) per side).
Rows: lemma degree floor d_min = 1/x - 1 (Caro-Wei on the correction graph); parking floor 1/(2a) (<= one width-1
helper per source); per-level counting (levels.py) = 0; design-side private-helper ceiling at degree d;
E8 per-correction heuristic (4142 helper mass / 5824 corrections per stage)."""
import sys, math
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # kappa_pair.py (expNB), transfer.py (expNE)
import kappa_pair as KP
from transfer import eta3

RHO_T = {'bit': 3512 / 1554, 'cx': 4142 / 1792}
OFF = {'bit': 0.1319 / eta3(8, 63 / 111, RHO_T['bit']), 'cx': 0.1350 / eta3(9, 64 / 112, RHO_T['cx'])}
EFF = 4142 / (112 * 52)          # E8 complex: helper width units per correction pair, per stage


def side(sd, h, x, rho, best=True):
    s3 = OFF[sd] * eta3(h, x, rho) * KP.X(3, h, x)
    return max(s3, KP.side_saving(sd, h, x, 'cal9')) if best else s3


def mult(h, x, rho, best=True):          # unit in complex orientation at h, bit twin h-1, same x, same rho both sides
    A = min(side('bit', h - 1, x, rho, best), side('cx', h, x, rho, best))
    return A / (1 + A) / 0.0017


def priv(d, a):                         # private helpers, target collects in-sources along its flag: widths 1..d (cap a)
    lo = math.floor(d); f = d - lo
    g = lambda k: sum(min(i, a) for i in range(1, k + 1))
    return ((1 - f) * g(lo) + f * g(lo + 1)) / (2 * a)


if __name__ == '__main__':
    a8 = 8
    print('E8 complex: rho_today = %.3f; parking floor 1/(2a) = %.4f; expSH relaxation LB 139/1792 = %.4f;'
          ' private-helper design at d=52: %.1f; EFF = %.3f width/correction' % (RHO_T['cx'], 1 / 16, 139 / 1792,
                                                                                   priv(52, a8), EFF))
    assert priv(52, 8) > RHO_T['cx'] > 139 / 1792 > 0      # bounds bracket today's word
    print('lemma minimum degree at E8 x=64/112: d_min = %.3f (actual 52 = every free pair)' % (112 / 64 - 1))
    xs = (0.30, 0.35, 0.40, 0.45, 0.50)
    for h in (5, 6):
        a = h - 1
        print('\nh_c = %d (a = %d)   x: ' % (h, a) + '  '.join('%5.2f' % x for x in xs))
        print('  d_min = 1/x-1           ' + '  '.join('%5.2f' % (1 / x - 1) for x in xs))
        print('  rho_LB (lemma+frames)   ' + '  '.join('%5.3f' % (1 / (2 * a)) for x in xs) + '   (<= parking floor; counting adds 0)')
        print('  rho private @ d_min     ' + '  '.join('%5.3f' % priv(1 / x - 1, a) for x in xs))
        for V in (2 ** (h - 1), 2 ** h, 2 ** (h + 1)):
            print('  rho heur dense V=%-4d    ' % V + '  '.join('%5.2f' % (EFF * (V / 2) / (2 * a)) for x in xs))
        for lab, rf in (('x at rho_LB', lambda x: 1 / (2 * a)), ('x at rho_priv(d_min)', lambda x: priv(1 / x - 1, a)),
                        ('x at rho=1.5', lambda x: 1.5), ('x at rho_today', lambda x: RHO_T['cx'])):
            print('  mult best %-20s' % lab + '  '.join('%5.2f' % mult(h, x, rf(x)) for x in xs)
                  + '   | s3 only ' + ' '.join('%5.2f' % mult(h, x, rf(x), False) for x in xs))
    print('\nrho at which s3 beats cal9 (h, x): ' + ', '.join(
        '(%d,%.2f): %s' % (h, x, next(('%.2f' % r for r in [i / 100 for i in range(400, 0, -1)]
                                      if side('cx', h, x, r, False) >= KP.side_saving('cx', h, x, 'cal9')), 'none'))
        for h in (5, 6) for x in (0.3, 0.45)))
    print('PASS_RHO_MIN')
