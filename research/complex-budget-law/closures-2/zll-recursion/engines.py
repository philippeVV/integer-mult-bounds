"""Pinned repo pricing engines, imported unchanged (slim copy of the scouting wrapper expM8/price.py).

complex: research/complex-local-word-ceiling/ceiling.py  (m=73, 896 families, 9 stages, idle cap 51, centre 64)
bit    : research/bit-revision3-levers-seed03/levers.py  (PR589: m=64, w=888/9, cst 63 = 9 rank-7 copies, cap 45)
The repository root is found from this file's location (research/complex-budget-law/closures-2/<lane>/);
set IMB_REPO to override.  kappa = a/(1+a) with a = min(complex cap, bit cap) (planning ceilings).
"""
import os
import sys
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(os.environ.get('IMB_REPO') or Path(__file__).resolve().parents[4])
sys.path.insert(0, str(ROOT / 'research' / 'complex-local-word-ceiling'))
sys.path.insert(0, str(ROOT / 'research' / 'bit-revision3-levers-seed03'))
import ceiling as Cx  # noqa: E402
import levers as Bt   # noqa: E402

CW = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82}                    # complex non-centre local profile per call (#1094)
BW = {1: 1323, 2: 1342, 3: 79, 4: 44, 5: 102, 6: 11, 7: 10}     # bit #589 local minus the 9 rank-7 copies
REF = Q('0.0017')
G = Q(1, 10 ** 9)


def cx(copies8, centre):
    loc = dict(CW)
    if copies8:
        loc[8] = copies8
    try:
        return Cx.ceiling(loc, 73, 896, 9, centre, grid=G)[1]
    except ValueError:
        return Q(0)


def bt(copies7, cst):
    loc = dict(BW)
    loc[7] = loc[7] + copies7
    H = Bt.relaxed(loc, 64, Q(888, 9), cst, 45)
    try:
        return Bt.root(H, 64, Q(888, 9), tol=G)[1]
    except Bt.Reject:
        return Q(0)


def kap(a):
    return a / (1 + a)
