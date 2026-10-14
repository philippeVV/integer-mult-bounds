"""Break-even copy charge for plane ports by decoder term rank, priced with #1240 ceiling.py (unpaid=2).
Local word = selected word with its 8 rank-8 centre copies replaced by T decoder terms of rank rho
(C = T*rho per call); rest of the paid local profile unchanged (first-order, as in #1587)."""
import sys
from fractions import Fraction as Q
sys.path.insert(0, './research/complex-local-word-ceiling')
import ceiling as C
WORD = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}
M, FAM, ST = 73, 896, 9
G = Q(1, 10**9)
base = C.ceiling(WORD, M, FAM, ST, 64, grid=G)[1]
print('baseline a_cap', float(base))
def cap(rho, T):
    loc = dict(WORD); del loc[8]
    loc[rho] = loc.get(rho, 0) + T
    try:
        return C.ceiling(loc, M, FAM, ST, rho * T, grid=G, unpaid=2)[1]
    except ValueError as e:
        return None
for rho in [1, 2, 3, 4, 6, 7, 8, 9]:
    lo, hi = 1, 400
    # largest T (fractional allowed? keep integer) with cap > base
    best = None
    while lo <= hi:
        mid = (lo + hi) // 2
        c = cap(rho, mid)
        if c is not None and c > base:
            best = mid; lo = mid + 1
        else:
            hi = mid - 1
    print(f'rho={rho}: max terms T={best}  max C={best*rho if best else None}  cap at T: {float(cap(rho,best)) if best else None}', flush=True)
