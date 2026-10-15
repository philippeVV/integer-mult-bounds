#!/usr/bin/env python3
"""Favorable price gate for the cheapest information-sufficient QQ16 final-record contract found by
modelq2.py (INIT64 + per-bank d_A on banks 0,1,2), against the current best claim (PR2474).

All numbers are ideal-packing ceilings from the pinned research/complex-local-word-ceiling/ceiling.py
(merged PR1240); a ceiling is an upper bound for the stated model, not an attained exponent.
Premise granted for free: the QQ16 controller removes the centre path (C = 0) and the local words
are the retained B/F profiles (5 F + 4 B calls per replica).

Retained arrays and lifetimes (calls during which each array must stay live, until the final decode):
  bank k first call stage: 0, 3, 5, 7.  INIT records of bank k are read just before that call.
  INIT64: 16 arrays per bank, lifetime 9, 6, 4, 2 calls          -> 16*21 = 336 array-calls
  d_A   : 8 arrays on banks 0,1,2, lifetime 8, 5, 3 calls          ->  8*16 = 128 array-calls
Conventions per array-call (n = array-calls / 9 per call):
  storage : centre += n (one paid unit per live array per call; extra mass packed as ideal idle
            children) -- the most favorable;
  centre  : n rank-8 children + centre += 8n (price3 / PR2470 convention: the selected word's
            own price of one centre query coordinate per call).
Optional known-construction surcharge (not a lower bound): each 8-coordinate initial query priced
as PR2468's native B2C read invoice (192 rank-one children), 8 queries per replica, plus one rank-one
line per retained d_A array (24).
Prepared by Philippe Van Velzen with Claude Code (Claude Opus 5.5). Apache-2.0.
"""
import json
import sys
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0, str(__import__('pathlib').Path(__import__('os').environ.get('REPO', __import__('pathlib').Path(__file__).resolve().parents[3])) / 'research/complex-local-word-ceiling'))
import ceiling as C  # noqa: E402

M, FAM, ST = 73, 896, 9
GRID = Q(1, 10 ** 12)
B = {1: 2436, 2: 418, 3: 432, 4: 142, 5: 82}
F = {1: 2420, 2: 418, 3: 432, 4: 142, 5: 82}
MIX = {r: Q(5 * F[r] + 4 * B[r], 9) for r in B}
SEL = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}
CLAIM = Q('0.001703676387993676')            # PR2474 conditional kappa (truncated)
BIT = Q('0.001733427039589844234749')        # PR589 bit supplier saving
INIT_AC = 16 * (9 + 6 + 4 + 2)
DA_AC = 8 * (8 + 5 + 3)


def cap(local, centre):
    try:
        return C.ceiling(local, M, FAM, ST, centre, grid=GRID)[1]
    except ValueError:
        return None


def kap(a):
    return None if a is None else a / (1 + a)


def scenario(ac, kind, reads=False):
    n = Q(ac, 9)
    p, c = dict(MIX), Q(0)
    if kind == 'storage':
        c += n
    elif kind == 'centre':
        p[8] = p.get(8, 0) + n
        c += 8 * n
    if reads:
        p[1] = p[1] + Q(8 * 192 + 24, 9)
    return cap(p, c)


def breakeven(kind, target_kappa):
    lo, hi = 0, 20000
    while hi - lo > 1:
        mid = (lo + hi) // 2
        a = scenario(mid, kind)
        if a is not None and kap(a) > target_kappa:
            lo = mid
        else:
            hi = mid
    return lo


def main():
    out = {}
    base = cap(MIX, 0)
    sel = cap(SEL, 64)
    out['free_controller_C0'] = str(kap(base))
    out['selected_word_C64'] = str(kap(sel))
    chk = cap({**MIX, 8: Q(8)}, 64)
    assert chk == scenario(72, 'centre')     # 8 carriers per call re-priced as the selected centre path
    print(f'control: QQ16 words + the selected centre path (8 rank-8 copies, C=64) kappa_cap {float(kap(chk)):.9f} '
          f'= centre convention at 72 array-calls')
    print(f'free QQ16 controller, C=0: kappa_cap {float(kap(base)):.6f}; selected word C=64: {float(kap(sel)):.9f}')
    print(f'claim PR2474 {float(CLAIM):.12f}; bit side {float(BIT):.12f}')
    ac = INIT_AC + DA_AC
    print(f'retained array-calls per replica: INIT64 {INIT_AC} + d_A {DA_AC} = {ac}')
    for kind in ('storage', 'centre'):
        for reads in (False, True):
            a = scenario(ac, kind, reads)
            k = kap(a)
            s = 'no positive saving' if k is None else f'kappa_cap {float(k):.9f} ({float(k / CLAIM - 1) * 100:+.2f}% vs claim)'
            print(f'  {kind:8s} reads={str(reads):5s}: {s}')
            out[f'{kind}_reads{reads}'] = None if k is None else str(k)
    for kind in ('storage', 'centre'):
        be = breakeven(kind, CLAIM)
        print(f'break-even array-calls per replica vs claim ({kind}): {be}')
        out[f'breakeven_{kind}'] = be
    Path(__file__).with_name('PRICE-GATE.json').write_text(json.dumps(out, indent=1))
    print('PASS_QQ16_FINAL_RECORD_PRICE_GATE')


if __name__ == '__main__':
    main()
