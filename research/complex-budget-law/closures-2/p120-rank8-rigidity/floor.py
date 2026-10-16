#!/usr/bin/env python3
"""Transfer-degree floor for a masked decoder with copy fraction gamma (exact, stdlib).

(i)  Caro-Wei on the transfer graph T: on a T-independent set Q, N[Q,Q] = I, so alpha(T) <= rank P+ <= n+, the number
     of positive-cost totals; hence mean degree d >= V/n+ - 1.
(ii) Hadamard capacity (expN2 L1): V <= sum_j (2^beta_j - 1) over totals of cost beta_j.
(iii) C = sum_j beta_j = gamma V.
LP over y_b = (# totals of cost b)/V, b = 1..12:  maximise sum y_b  s.t.  sum b y_b = gamma,  sum (2^b - 1) y_b >= 1.
Two constraints, so an optimum sits at a vertex with <= 2 nonzero y_b; we enumerate every such vertex exactly.
Also prints the 'full-cost' floor h/gamma - 1 (every total costs h).
Run: python3 -I -B floor.py   (prints PASS_FLOOR)
"""
from fractions import Fraction as Q
from itertools import combinations

BS = range(1, 13)


def lp(g):
    best = None
    cands = []
    for b in BS:                                       # one nonzero: y_b = g/b, needs (2^b-1) g/b >= 1
        y = g / b
        if (2 ** b - 1) * y >= 1:
            cands.append(({b: y}, y))
    for b1, b2 in combinations(BS, 2):                 # two nonzero, both constraints tight
        # b1 y1 + b2 y2 = g ; (2^b1-1) y1 + (2^b2-1) y2 = 1
        a11, a12, a21, a22 = b1, b2, 2 ** b1 - 1, 2 ** b2 - 1
        det = a11 * a22 - a12 * a21
        if det == 0:
            continue
        y1 = (g * a22 - a12 * 1) / det
        y2 = (a11 * 1 - a21 * g) / det
        if y1 >= 0 and y2 >= 0:
            cands.append(({b1: y1, b2: y2}, y1 + y2))
    for mix, val in cands:
        if best is None or val > best[1]:
            best = (mix, val)
    return best


if __name__ == '__main__':
    print('gamma   max n+/V   d_floor(CW+Hadamard)  optimal cost mix (y_b)          | full-cost floor h/gamma-1: h=9  h=12')
    ok = True
    for g in (Q(3, 10), Q(2, 5), Q(1, 2), Q(4, 7), Q(3, 5), Q(7, 10), Q(4, 5), Q(8, 9)):
        mix, val = lp(g)
        d = 1 / val - 1
        print(f'{float(g):.3f}   {float(val):.4f}     {float(d):6.3f}               '
              f'{ {b: round(float(y), 5) for b, y in mix.items()} }'.ljust(70) + f'|  {float(9 / g - 1):5.2f}  {float(12 / g - 1):5.2f}')
        # the CW+Hadamard floor never forces dense transfers: it stays below 1/gamma (cost-1 totals dominate)
        ok &= d < 1 / g and 9 / g - 1 > d
    g = Q(4, 7)
    print(f'E8 (gamma = 4/7): floor {float(1 / lp(g)[1] - 1):.3f} (~1/gamma - 1 = 0.75); full-cost floor 14.75 vs realised 52-56')
    print('PASS_FLOOR' if ok else 'FAIL_FLOOR')
