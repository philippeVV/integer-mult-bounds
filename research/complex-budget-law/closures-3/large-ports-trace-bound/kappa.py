"""expNM: kappa planning ceilings (expM1 framework_ceiling, expL4 large_port form, data-floor regime, best odd s)
for each family's (h, u) at gamma = best exact (1), LM-dense (minrank*h/(uV)), and the proven floor. Stdlib."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'framework-ceiling'))
from framework_ceiling import A_cap, kap, lbar_of
def lp(h, u, g):
    best = (0.0, None)
    for s in (3, 5, 7, 9):
        m = u + s * (h - 1); nu = u * (1 - s / (s - 1) * g)
        if nu <= 0: continue
        kk = kap(A_cap(m, nu, lbar_of(m, min(1.0, 2 / (s - 1)), h - 1)))
        if kk > best[0]: best = (kk, s)
    return best
rows = [  # name, h, u, V, LM minrank (numerical), alpha floor
 ('PG(2,4) dot', 6, 2, 16, 7, 6), ('PG(3,4) herm', 8, 2, 40, 10, 7), ('PG(3,4) dot', 8, 2, 64, 16, 10),
 ('AG(3,4) dot +1', 8, 2, 48, 16, 10), ('AG(2,8) dot +1', 9, 3, 56, 24, 16),
 ('PG(2,16) herm', 12, 4, 208, 90, 44), ('PG(2,16) dot (frames, LM not run)', 12, 4, 256, 92, 52), ('PG(3,8) dot (frames)', 12, 3, 512, 143, 40)]
print(f"{'family':16s} {'h':>3s} {'u':>2s} {'k':>5s} {'V':>4s} | gamma: exact  LMdense  floor | ceiling x over 0.0017 at each (s)")
for nm, h, u, V, K, al in rows:
    g = {'exact': 1.0, 'LMdense': K * h / (u * V), 'floor': al / V}
    cells = []
    for key in ('exact', 'LMdense', 'floor'):
        kk, s = lp(h, u, g[key]); cells.append(f"{key} {g[key]:.3f}: " + (f"x{kk / 0.0017:5.1f} (s{s})" if s else 'nu<=0 (no gain)'))
    print(f"{nm:16s} {h:3d} {u:2d} {u / (h - 1):5.3f} {V:4d} | " + ' | '.join(cells))
assert all(lp(h, u, 1.0)[1] is None for _, h, u, *_ in rows)   # gamma = 1 gives nu <= 0 everywhere
print('x10-calibrated needs ceiling ~x26: gamma <= 0.21 (k=1/4) or 0.31-0.34 (k=3/8) [expL4 s9]; stop test: best exact > 0.705 everywhere')
print('PASS_EXPNM_KAPPA')
