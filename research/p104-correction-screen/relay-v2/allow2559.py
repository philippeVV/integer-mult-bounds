import json, os, sys
from fractions import Fraction as Q
ROOT = sys.argv[1]; sys.path.insert(0, os.path.join(ROOT, 'research', 'complex-local-word-ceiling'))
import ceiling as CE
oc = json.load(open(sys.argv[2])); pr = json.load(open(sys.argv[3]))
k0 = Q(oc['kappa_rational']) if '/' in str(oc['kappa_rational']) else None
if k0 is None:
    s = oc['kappa']; ip, fp = s.split('.'); k0 = Q(int(ip + fp), 10**len(fp))
bitc = Q(pr['prices']['native']['outer']['bit_c'])
k = Q(-((-k0.numerator * 10**30) // k0.denominator), 10**30)
print('kappa', float(k), 'bitcap', float(bitc))
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
    print(p, maxN(p, k), maxN(p, bitc), flush=True)
