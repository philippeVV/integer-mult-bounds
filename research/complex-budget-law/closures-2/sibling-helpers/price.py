#!/usr/bin/env python3
"""expNU / price.py -- idea 8 (sibling calls share helper registers), priced on BOTH suppliers.

Engines (pinned repo code, imported unchanged, same calls as expM8/price.py):
  complex: research/complex-local-word-ceiling/ceiling.py  (m=73, 896 families, 9 stages, idle cap 51, centre 64)
  bit    : research/bit-revision3-levers-seed03/levers.py  (PR589: m=64, w=888/9, cst 63, idle cap 45)
Per-call local profiles are split into data + centre copies + helper climbs:
  complex #1094: helper climbs {1:1675,2:272,3:365,4:54,5:40} (480 helpers, mass 3730; wordmodel replay, expSK)
  bit #589     : helper runs   {1:994,2:1041,3:70,4:44,5:10}  (439 helpers, mass 3512; expNE bit_inputs / #555 split)
Sharing model: a fraction sigma of the per-call helper climb children is served jointly by k sibling calls, so
each call pays (1-sigma) + sigma/k of today's helper children; the freed coordinates go back to ideal idle
(coverage law, engine-internal).  Junk handling is charged as an optional per-shared-helper rank-1 buffer.
kappa = a/(1+a), a = min(complex cap, bit cap); complex-only rows are capped by the realised bit endpoint.
"""
import os
import sys
from fractions import Fraction as Q
from pathlib import Path
ROOT = Path(os.environ.get('IMB_REPO') or Path(__file__).resolve().parents[4])   # repository root
sys.path.insert(0, str(ROOT / 'research' / 'complex-local-word-ceiling'))
sys.path.insert(0, str(ROOT / 'research' / 'bit-revision3-levers-seed03'))
import ceiling as Cx  # noqa: E402
import levers as Bt   # noqa: E402

G = Q(1, 10 ** 9)
REF, OURS, BIT_REAL = Q('0.0017'), Q('0.0017148565'), Q('0.0017334270395898442')
CX_H = {1: 1675, 2: 272, 3: 365, 4: 54, 5: 40}
CX_REST = {1: 737, 2: 146, 3: 67, 4: 88, 5: 42, 8: 8}           # data climbs + 8 rank-8 centre copies
BT_H = {1: 994, 2: 1041, 3: 70, 4: 44, 5: 10}
BT_LOC = {1: 1323, 2: 1342, 3: 79, 4: 44, 5: 102, 6: 11, 7: 10 + 9}   # #589 local incl. 9 rank-7 copies
BT_REST = {r: BT_LOC[r] - BT_H.get(r, 0) for r in BT_LOC}
NH = {'cx': 480, 'bit': 439}


def mass(d):
    return sum(r * n for r, n in d.items())


def prof(rest, H, scale, extra1=0):
    p = {r: Q(n) for r, n in rest.items()}
    for r, n in H.items():
        p[r] = p.get(r, 0) + Q(n) * scale
    p[1] = p.get(1, 0) + extra1
    return {r: n for r, n in p.items() if n > 0}


def cx(p, centre=64):
    try:
        return Cx.ceiling(p, 73, 896, 9, centre, grid=G)[1]
    except ValueError:
        return Q(0)


def bt(p, cst=63):
    try:
        return Bt.root(Bt.relaxed(p, 64, Q(888, 9), cst, 45), 64, Q(888, 9), tol=G)[1]
    except Bt.Reject:
        return Q(0)


def caps(sigma, k, junk=0):
    """junk: rank-1 buffer children per shared helper per call (expSL (b) convention)."""
    s = (1 - sigma) + sigma / Q(k) if k else (1 - sigma)
    jc = junk * sigma * NH['cx'] / (k or 1)
    jb = junk * sigma * NH['bit'] / (k or 1)
    return cx(prof(CX_REST, CX_H, s, jc)), bt(prof(BT_REST, BT_H, s, jb))


def kap(a):
    return a / (1 + a)


def main():
    assert mass(CX_REST) + mass(CX_H) == 5586 and mass(BT_LOC) == 5129
    assert mass(BT_REST) == 1554 + 63 and mass(BT_H) == 3512
    c0, b0 = caps(0, 1)
    print(f'baseline caps: complex {float(c0):.7f} bit {float(b0):.7f}  (expM8: 0.0016901 / 0.0019559)')
    print(f'helper share of local mass: complex {mass(CX_H)}/5586 = {mass(CX_H)/5586:.3f};'
          f' bit {mass(BT_H)}/5129 = {mass(BT_H)/5129:.3f}')
    print('\nsigma = shared fraction of helper children, k = sibling calls per shared helper, junk = rank-1 buffers per shared helper')
    print(f'{"sigma":>5} {"k":>3} {"junk":>4} {"cx cap":>10} {"x":>6} {"bit cap":>10} {"x":>6} {"kappa(both)":>11} {"/0.0017":>7} {"cx-only kappa":>13} {"/0.0017":>7}')
    rows = [(1, 2, 0), (1, 2, 1), (1, 4, 0), (1, 0, 0), (Q(1, 2), 2, 0), (Q(1, 4), 2, 0), (Q(1, 10), 2, 0),
            (Q(1, 2), 0, 0), (Q(1, 4), 0, 0)]
    for sigma, k, junk in rows:
        a, b = caps(sigma, k, junk)
        both = kap(min(a, b))
        cxo = kap(min(a, BIT_REAL))
        print(f'{float(sigma):5.2f} {k if k else "inf":>3} {junk:4d} {float(a):10.7f} {float(a/c0):6.3f} {float(b):10.7f} '
              f'{float(b/b0):6.3f} {float(both):11.7f} {float(both/REF):7.3f} {float(cxo):13.7f} {float(cxo/REF):7.3f}')
    # break-even: how much shared fraction gives +5% / +10% on the BOTH-sides kappa at k=2
    print('\nk=2, junk=0: smallest sigma (grid 1/100) reaching a both-sides kappa gain over today\'s caps:')
    base = kap(min(c0, b0))
    for tgt in (Q(105, 100), Q(110, 100), Q(125, 100)):
        hit = next((s for s in range(0, 101) if kap(min(*caps(Q(s, 100), 2))) >= tgt * base), None)
        print(f'  x{float(tgt):.2f}: sigma >= {hit/100 if hit is not None else "none"}')
    print('PASS_EXPNU_PRICE')


if __name__ == '__main__':
    main()
