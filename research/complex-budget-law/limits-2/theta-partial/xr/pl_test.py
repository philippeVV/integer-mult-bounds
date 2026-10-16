"""Do the point-flat or line-flat gauges (Theorem P and its mirror) alone certify Theta3 <= 1/8?"""
import sys, numpy as np
from model import *
from xk_sparse import flat_gauge
import xg_core as Dc
name = sys.argv[1]; V, nP, E, g = load_any(name); S = Struct(V, nP, E); rng = np.random.default_rng(0)
codes = [('rand%d' % k, Dc.rand_codes(V, rng)) for k in range(6)]
for f in sys.argv[2:]: codes.append((f.split('/')[-1], [orth(np.asarray(q, complex)) for q in np.load(f, allow_pickle=True)]))
for lab, Q in codes:
    xp = flat_gauge(S, Q, 0, 3000)[0]; xl = flat_gauge(S, Q, 1, 3000)[0]
    tp = theta(S, Q, 3, xp); tl = theta(S, Q, 3, xl)
    print(f'{lab}: point-flat {tp:.7f}  line-flat {tl:.7f}  min {min(tp, tl):.7f} {"<= 1/8" if min(tp, tl) <= 0.125 else "> 1/8"}', flush=True)
