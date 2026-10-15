"""Partial-upgrade break-even under #1240 ceiling.py: unpaid = 1 + |S|/112, centre = 64 + dC, extra dC terms of rank rho."""
import sys
from fractions import Fraction as Q
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../complex-local-word-ceiling'))
import ceiling as C
WORD = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}
M, FAM, ST = 73, 896, 9; G = Q(1, 10**10)
base = C.ceiling(WORD, M, FAM, ST, 64, grid=G)[1]
print('baseline cap', float(base))
for s in (8, 16):
    for dC in range(0, 2 * s + 1):
        for rho in (1, 2):
            if dC % rho: continue
            loc = dict(WORD); loc[rho] = loc.get(rho, 0) + dC // rho
            try: cap = C.ceiling(loc, M, FAM, ST, 64 + dC, grid=G, unpaid=1 + Q(s, 112))[1]
            except ValueError: cap = None
            print(f'|S|={s} dC={dC} extra terms rank {rho}: cap {float(cap):.10f} {"WIN" if cap > base else "lose"} ratio {float(cap/base):.5f}')
