"""EXACT universal floor alpha(Gamma_E) >= 14 for EVERY plane assignment (hence C >= 15 for every decoder).
W = j^perp (even-weight vectors); dot is alternating there. For a Lagrangian L (4-dim isotropic) put S_L = {j+w : 0 != w in L}:
odd vectors with pairwise dot 1. Ports whose plane meets S_L are pairwise non-orthogonal (x.y = 1, or x = y odd).
No plane contains j (r.j = r.r = 1 != 0); a plane holds at most one vector of S_L (its two odd vectors are orthogonal).
Every nonzero w lies in exactly 135 of the 2295 Lagrangians  =>  average #ports hitting S_L = 224*135/2295 = 13.18."""
import json, sys, os
from collections import Counter
from common import *
J = 511
EVEN = [x for x in range(1, 512) if bin(x).count('1') % 2 == 0]
def span(B):
    s = {0}
    for b in B: s |= {x ^ b for x in s}
    return frozenset(s)
lag = set()
def grow(B, sp):
    if len(B) == 4: lag.add(sp); return
    for x in EVEN:
        if x > max(B, default=0) and x not in sp and all(dot(x, b) == 0 for b in B):
            grow(B + [x], span(B + [x]))
grow([], frozenset({0}))
cnt = Counter(w for L in lag for w in L if w)
print('Lagrangians', len(lag), 'per nonzero even vector', set(cnt.values()), 'vectors covered', len(cnt))
assert len(lag) == 2295 and set(cnt.values()) == {135} and len(cnt) == 255
print('average ports hitting S_L for ANY r:', 224 * 135, '/', 2295, '=', round(224 * 135 / 2295, 3), '=> alpha(Gamma_E) >= 14 for all r')
for f in sys.argv[1:]:
    d = json.load(open(f)); runs = d if isinstance(d, list) else [d]
    for run in runs[:1]:
        r = run['r']; E = planes(r); best = 0
        for L in lag:
            S = {J ^ w for w in L if w}
            I = [a for a in range(n) if P[a] in S or r[a] in S]
            assert all(not orth_pl(E[a], E[b]) for a in I for b in I if a < b)
            best = max(best, len(I))
        print(os.path.basename(f), 'best Lagrangian independent set', best)
