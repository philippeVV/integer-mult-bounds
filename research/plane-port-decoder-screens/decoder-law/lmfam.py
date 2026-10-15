"""LM on a fixed P120 family of symmetric hyperplane-section supports (h list as JSON), many restarts."""
import sys, json
import numpy as np
from p120lm import *
hs = json.loads(sys.argv[1]); tries = int(sys.argv[2])
sup = [(HYP[h], HYP[h]) for h in hs]
print('cost', cost(sup), 'boolsat sym', boolsat(sup, sym=True), 'asym', boolsat(sup), flush=True)
for sym in (True, False):
    best = 1e9
    for s in range(tries):
        r = lm(FULL, *masks(sup), seed=100 + s, tries=1, maxfev=600, sym=sym)
        best = min(best, r[0])
        if r[0] < 1e-9: print('EXACT-ish', sym, s); break
    print('sym' if sym else 'asym', 'best residual', best, 'squared', round(best ** 2, 3), flush=True)
