"""Pair scoring at h_c = 5 (bit twin h_b = 4 by the expNB rule, x_b = x_c): expNB 'cal' and expNE 'best' rows."""
import os, sys, io, contextlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # pair_s3.py, transfer.py (expNE), kappa_pair.py (expNB)
with contextlib.redirect_stdout(io.StringIO()):
    import pair_s3 as P
import kappa_pair as KP
for x in [0.3125, 0.4027, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]:
    cal = KP.kappa_pair(4, x, x, 1, 'cal')[0]; best = P.pair(5, x, 'best')
    print('h_c=5 x=%.4f  cal %.5f (x%.2f)  best %.5f (x%.2f)' % (x, cal, cal / 0.0017, best, best / 0.0017))
# break-even x (pair = 0.0017) on each row
def be(f):
    lo, hi = 0.3, 0.99
    for _ in range(40):
        mid = (lo + hi) / 2
        if f(mid) > 0.0017: lo = mid
        else: hi = mid
    return lo
print('break-even x (x1 over 0.0017): cal %.4f  best %.4f' % (be(lambda x: KP.kappa_pair(4, x, x, 1, 'cal')[0]), be(lambda x: P.pair(5, x, 'best'))))
print('x for x2: cal %.4f best %.4f' % (be(lambda x: KP.kappa_pair(4, x, x, 1, 'cal')[0] / 2), be(lambda x: P.pair(5, x, 'best') / 2)))
assert abs(P.pair(5, 0.4027, 'best') / 0.0017 - 3.03) < 0.01 and be(lambda x: P.pair(5, x, 'best')) < 0.73
print('PASS_EXPNF_SCORE')
