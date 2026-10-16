#!/usr/bin/env python3
"""expNE transfer: eta3(h, x; rho) for a small unit priced with E8's per-port run statistics at s = 3.
Per port per stage: data mass 2a in today's absolute run-size mix (bit #555 x+y runs, ranks capped at a, mass kept),
helper mass rho * 2a in today's helper run mix (#555), copies x at cap h.  No idle (K = 1), helpers on exactly tiled
banks (D unchanged).  eta3 = root / X(3, h, x).  Planning model only: assumes the small unit's helper mass per data
mass (rho) and fragmentation equal E8's (rho_today = 3512/1554 = 2.26 bit, 4142/1792 = 2.31 complex)."""
import math
from collections import Counter

DATA = {1: 327, 2: 301, 3: 10, 4: 3, 5: 92, 6: 10, 7: 9}          # #555 x+y per stage, mass 1554 = 2*111*7
HELP = {1: 994, 2: 1041, 3: 70, 4: 44, 5: 10}                    # #555 helper runs, mass 3512


def root(prof, m, D):
    f = lambda q: sum(z * r * math.expm1(q * math.log(m / r)) for r, z in prof.items() if r < m) - D
    lo, hi = 0.0, 0.5
    for _ in range(100):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return hi


def capped(mix, a, total):
    c = Counter()
    for r, z in mix.items():
        c[min(r, a)] += r * z / min(r, a)
    s = sum(r * z for r, z in c.items())
    return {r: z * total / s for r, z in c.items()}


def eta3(h, x, rho, frag=True):
    a, m = h - 1, 3 * h - 2
    D = 2 - 3 * x
    X3 = root({a: 6, h: 3 * x / h}, m, D)
    prof = Counter()
    for r, z in (capped(DATA, a, 2 * a) if frag else {a: 2.0}).items():
        prof[r] += 3 * z
    for r, z in capped(HELP, a, rho * 2 * a).items():
        prof[r] += 3 * z
    prof[h] += 3 * x / h
    return root(prof, m, D) / X3


if __name__ == '__main__':
    print('eta3 = realised/X(3) (rows h, cols rho = helper mass / data mass; x = 0.30 / 0.45 / 0.5676)')
    for x in (0.30, 0.45, 63 / 111):
        print('x = %.4f' % x)
        print('  h   ' + ''.join('rho=%-6.2f' % r for r in (0, 0.5, 1.09, 1.5, 2.26)) + ' | data at cap a, rho=2.26')
        for h in (5, 6, 7, 8, 9):
            print('  %d   ' % h + ''.join('%-10.3f' % eta3(h, x, r) for r in (0, 0.5, 1.09, 1.5, 2.26))
                  + ' | %.3f' % eta3(h, x, 2.26, frag=False))
    e = eta3(8, 63 / 111, 3512 / 1554)
    print('check: E8 bit row h=8, x=63/111, rho=2.26 -> %.4f (price_s3 exact-tiled 0.1503)' % e)
