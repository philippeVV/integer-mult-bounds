"""expNV / price_flat.py -- price the flat tower (n > 4 pairs, m fixed) with the repo ceiling engines, both suppliers.
Engines: engines.py (pinned research/complex-local-word-ceiling/ceiling.py and bit-revision3-levers-seed03/levers.py),
formula exactly as expNL/price_n.py: complex ceiling(local, 73, 224n, units, 64 + extra), bit relaxed/root at
m=64, w=222n/units, cst 63 + extra.  kappa = a/(1+a), a = min(complex, bit).

Rows
 D  same-family stage sharing (2n+1 calls, m fixed): proved repair = descents.  Along the call sequence the frame
    dimension must rise 8 per call (complex; 7 bit) but only 9 blocks exist, so every port register chain descends
    >= 8*2(n-4) (bit 7*2(n-4)) per replica, and every descended unit is re-climbed (coverage: net climb fixed).
    Extra mass per replica >= 2*descent = 2*16(n-4)*112 complex, 2*14(n-4)*111 bit.
      D-fav : all extra mass priced as ideal idle (rank R children) -- the most favourable admissible charge;
      D-real: descents as rank-8 (bit rank-7) children, re-climb as ideal idle.
 P  different-family sharing (parallel calls, no descent): impossible at 2n+1 units for n>=5 (levels.py theorem);
    with e extra units the counting bound gives e >= max(2, ceil(n/5)) (README 2c).  Ceiling only -- existence open.
Run: python3 -B price_flat.py
"""
import sys, os, math
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from fractions import Fraction as Q
import engines as P
Cx, Bt = P.Cx, P.Bt
def cx(fam, units, extra=Q(0), add8=Q(0), m=73):
    loc = dict(P.CW); loc[8] = 8 + add8
    try: return Cx.ceiling(loc, m, fam, units, 64 + extra, grid=P.G)[1]
    except (ValueError, AssertionError): return Q(0)
def bt(fam, units, extra=Q(0), add7=Q(0), m=64):
    loc = dict(P.BW); loc[7] = loc[7] + 9 + add7
    try:
        H = Bt.relaxed(loc, m, Q(fam, units), 63 + extra, Cx.max_idle_rank(m))
        return Bt.root(H, m, Q(fam, units), tol=P.G)[1]
    except Exception: return Q(0)
c0, b0 = cx(896, 9), bt(888, 9)
assert (c0, b0) == (P.cx(8, 64), P.bt(9, 63))
k0 = P.kap(min(c0, b0))
print(f'baseline n=4 selected: complex {float(c0):.7f} bit {float(b0):.7f} kappa {float(k0):.6f} (x{float(k0/P.REF):.3f} over 0.0017)')
def row(tag, n, units, ac, ab):
    a = min(ac, ab); k = P.kap(a)
    print(f'{tag:7} n={n:5} units={units:5}: complex {float(ac):.7f} (x{float(ac/c0):.3f})  bit {float(ab):.7f} '
          f'(x{float(ab/b0):.3f})  kappa {float(k):.6f}  x{float(k/P.REF):.3f} over 0.0017  x{float(k/k0):.3f} over n=4')
    return k
best = {}
print('\nD: same-family stage sharing, 2n+1 calls, m=73/64 fixed, paid descents (stop test fired: rank 8 per port)')
for n in (5, 6, 8, 16, 64, 1024):
    u = 2 * n + 1
    exc = Q(2 * 16 * (n - 4) * 112, u); exb = Q(2 * 14 * (n - 4) * 111, u)
    k1 = row('D-fav', n, u, cx(224 * n, u, exc), bt(222 * n, u, exb))
    k2 = row('D-real', n, u, cx(224 * n, u, exc, Q(2 * (n - 4) * 112, u)), bt(222 * n, u, exb, Q(2 * (n - 4) * 111, u)))
    best['D'] = max(best.get('D', 0), k1)
    w = Q(224 * n, u); wb = Q(222 * n, u)
    print(f'        coverage deficit test: positive saving needs extra < unpaid*w - centre: complex {float(exc):.1f} vs '
          f'{float(w - 64):.1f}, bit {float(exb):.1f} vs {float(wb - 63):.1f} -> engine rejects (a <= 0)')
    assert exc > w - 64 and exb > wb - 63
# break-even: largest extra mass per call (ideal idle) that keeps the n=5 flat tower at the n=4 baseline
lo, hi = Q(0), Q(48)
for _ in range(40):
    mid = (lo + hi) / 2
    if min(cx(1120, 11, mid), bt(1110, 11, mid)) >= min(c0, b0): lo = mid
    else: hi = mid
print(f'break-even at n=5: total repair mass <= {float(lo):.2f} per call (= {float(lo*11/112):.3f} rank-units per port per '
      f'replica) vs proved >= {3584/11:.1f} per call (16 per port): short by x{float(Q(3584,11)/lo):.0f}')
print('\nP: parallel calls on different families, no descent, e extra units (ceiling; existence not shown)')
for n in (5, 6, 10, 15, 20, 100, 1000):
    e = max(2, math.ceil(Q(n, 5))); u = 2 * n + e
    k = row('P-ceil', n, u, cx(224 * n, u), bt(222 * n, u)); best['P'] = max(best.get('P', 0), k)
print(f'\nbest D (favourable) kappa {float(best["D"]):.6f} = x{float(best["D"]/P.REF):.3f} over 0.0017')
print(f'best P ceiling   kappa {float(best["P"]):.6f} = x{float(best["P"]/P.REF):.3f} over 0.0017 '
      f'(complex-only cap with realised bit endpoint: {float(P.kap(min(Q("0.0017334"), cx(1120, 11)))):.6f})')
assert best['D'] < k0
print('PASS_EXPNV_PRICE')
