import json, os, sys
from fractions import Fraction as Q
sys.path.insert(0, '/home/pvv/repo/imb-work10/research/complex-local-word-ceiling')
import ceiling as CE
s = json.load(open('pr2559/OUTER-CLAIM.json'))['kappa']; ip, fp = s.split('.'); k0 = Q(int(ip + fp), 10**len(fp))
K = Q(-((-k0.numerator * 10**30) // k0.denominator), 10**30)
chan = {1: 1074, 2: 58, 3: 86, 4: 104, 5: 48, 6: 24, 7: 8, 8: 8}; m, fam, st, cen = 73, 832, 9, 56
def passes(prof, a=K):
    try: return CE.relax(prof, m, fam, st, cen, a)[1] < 1
    except ValueError: return False
def maxN(per, hi=2000):
    lo = 0
    while hi - lo > 1:
        mid = (lo + hi)//2
        o = dict(chan)
        for r, c in per.items(): o[r] = o.get(r, 0) + c*mid
        lo, hi = (mid, hi) if passes(o) else (lo, mid)
    return lo
if __name__ == '__main__':
    for p in sys.argv[1:]:
        steps = [int(x) for x in p.split(',')]; per = {}
        for r in steps: per[r] = per.get(r, 0) + 1
        print(p, maxN(per), flush=True)
