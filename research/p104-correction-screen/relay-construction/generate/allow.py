import json, os, sys
from fractions import Fraction as Q
ROOT = sys.argv[1]; sys.path.insert(0, os.path.join(ROOT, 'research', 'complex-local-word-ceiling'))
import ceiling as CE
led = json.load(open(sys.argv[2]))['rows']['struct_boundary_radix_union']
k0 = Q(led['kappa']['rational']); gap0 = Q(led['strict_checks']['actual_bit_above_effective']['rational'])
k = Q(-((-k0.numerator * 10**30) // k0.denominator), 10**30)
gap = Q((gap0.numerator * 10**30) // gap0.denominator, 10**30)
assert k >= k0 and k - k0 < Q(1, 10**29)
print('kappa', float(k), 'bitcap', float(k + gap))
chan = {1: 1074, 2: 58, 3: 86, 4: 104, 5: 48, 6: 24, 7: 8, 8: 8}; m, fam, st, cen = 73, 832, 9, 56
def add(b, e, kk):
    o = dict(b)
    for r, c in e.items(): o[r] = o.get(r, 0) + c*kk
    return o
def passes(prof, a):
    try: return CE.relax(prof, m, fam, st, cen, a)[1] < 1
    except ValueError: return False
def maxN(per, a, hi=2000):
    lo = 0
    while hi - lo > 1:
        mid = (lo + hi)//2; lo, hi = (mid, hi) if passes(add(chan, per, mid), a) else (lo, mid)
    return lo
for p in ({1: 2, 7: 1}, {1: 3, 6: 1}, {1: 4, 5: 1}, {1: 9}):
    print(p, maxN(p, k), maxN(p, k + gap), flush=True)
