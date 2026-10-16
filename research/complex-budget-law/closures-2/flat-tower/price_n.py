"""Planning ceilings for the only lever left by expNL: more exchange pairs n per replica.
Hosted typed model: n pairs need 2n+1 centre calls (proved lower bound, attained by the retained chain).
Families per call w = 2n*ports/(2n+1): complex 224n/(2n+1) (n=4: 896/9), bit 222n/(2n+1) (n=4: 888/9).
Variant F (favourable): m fixed (73 / 64), selected per-call local profile unchanged.
Variant T (tower): each call is one fresh-block stage, so m = 1+8(2n+1) complex, 1+7(2n+1) bit; idle cap from the
  halving guard; local profile unchanged.  Engines: engines.py (pinned ceiling.py / levers.py), unchanged.
Run: python3 -B price_n.py"""
import sys, os
here=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from fractions import Fraction as Q
import engines as P
Cx, Bt = P.Cx, P.Bt
def cx(n, m=73):
    loc=dict(P.CW); loc[8]=8
    try: return Cx.ceiling(loc, m, 224*n, 2*n+1, 64, grid=P.G)[1]
    except (ValueError, AssertionError) as e: return None
def bt(n, m=64):
    loc=dict(P.BW); loc[7]=loc[7]+9
    R=Cx.max_idle_rank(m)
    try:
        H=Bt.relaxed(loc, m, Q(222*n, 2*n+1), 63, R); return Bt.root(H, m, Q(222*n, 2*n+1), tol=P.G)[1]
    except Exception as e: return None
c0, b0 = P.cx(8,64), P.bt(9,63)
assert cx(4)==c0 and bt(4)==b0, (cx(4), c0, bt(4), b0)
print(f'baseline n=4 (selected template): complex {float(c0):.7f}  bit {float(b0):.7f}')
rows=[]
for n in (2,3,4,5,6,8,16,64,1024):
    for var in (('F','T') if n<=64 else ('F',)):
        mc = 73 if var=='F' else 1+8*(2*n+1); mb = 64 if var=='F' else 1+7*(2*n+1)
        ac, ab = cx(n, mc), bt(n, mb)
        if ac is None or ab is None:
            print(f'n={n:5d} {var}: engine rejects (m={mc}/{mb})'); continue
        a=min(ac,ab); k=P.kap(a)
        print(f'n={n:5d} {var}: w_cx {float(Q(224*n,2*n+1)):7.3f} m {mc:5d}/{mb:5d}  complex {float(ac):.7f} (x{float(ac/c0):.3f})  '
              f'bit {float(ab):.7f} (x{float(ab/b0):.3f})  kappa {float(k):.6f} x{float(k/P.REF):.3f} over 0.0017; '
              f'complex-only {float(P.kap(min(ac,Q("0.0017334")))):.6f}')
# the unhostable 2n-call ideal at n=4 (expND witness scaled): w = 896/8 = 112
loc=dict(P.CW); loc[8]=8
a8=Cx.ceiling(loc,73,896,8,64,grid=P.G)[1]
locb=dict(P.BW); locb[7]+=9
b8=Bt.root(Bt.relaxed(locb,64,Q(888,8),63,45),64,Q(888,8),tol=P.G)[1]
print(f'2n calls at n=4 (what the 4-centre witness would give if hostable): complex {float(a8):.7f} (x{float(a8/c0):.3f}) '
      f'bit {float(b8):.7f} (x{float(b8/b0):.3f}) kappa {float(P.kap(min(a8,b8))):.6f} x{float(P.kap(min(a8,b8))/P.REF):.3f}')
print('PASS_EXPNL_PRICE')
