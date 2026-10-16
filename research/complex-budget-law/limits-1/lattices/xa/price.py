"""expXA price.py: price gamma(h) = 6/(h-1) (LZD family, odd n = h) against the E8 control in the same formulas.

(1) expNM kappa.lp: expM1 framework ceiling, data-floor regime, best odd s in {3,5,7,9}, m = u + s(h-1).
(2) selected-word regime (expM1 section A): s = 9 template, local word share phi_sel / r_loc_sel at m.
(3) expQH peel: mu = 0 -> a = ln(m/(m-nu))/ln(m/(m-u)) ; mu = 1 -> exact root of {1:1, B: 1+(u-nu)/B}, B = m-u-1.
Multipliers are ceiling / 0.0017 and, more honestly, ratio to the E8 P112 control (h=9, gamma=4/7) in the same row.
"""
import sys, os, io, contextlib, math
HERE = os.path.dirname(os.path.abspath(__file__))
CBL = os.path.join(HERE, '..', '..', '..')          # research/complex-budget-law (in-branch copies of the scouting libraries)
sys.path.insert(0, os.path.join(CBL, 'closures-3', 'large-ports-trace-bound'))   # kappa.py (expNM)
sys.path.insert(0, os.path.join(CBL, 'framework-ceiling'))                       # framework_ceiling.py (expM1)
sys.path.insert(0, os.path.join(CBL, 'closures-3', 'peel-dead'))                 # hyp_drop.py (expQH)
with contextlib.redirect_stdout(io.StringIO()):
    from kappa import lp
from framework_ceiling import A_cap, kap, lbar_of
from hyp_drop import root

PHI_SEL, RLOC_SEL = None, None
word = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}
loc = sum(r * n for r, n in word.items())
PHI_SEL = loc / (72 * 896 / 9 + 64)
RLOC_SEL = 73 * math.exp(-sum(r * n * math.log(73 / r) for r, n in word.items()) / loc)


def sel(h, g, s=9, u=1):
    m = u + s * (h - 1)
    nu = u - s / (s - 1) * g * u
    return kap(A_cap(m, nu, lbar_of(m, PHI_SEL, RLOC_SEL)))


def peel(h, g, s=9, u=1, mu=0):
    m = u + s * (h - 1)
    nu = u * (1 - s / (s - 1) * g)
    if nu <= 0:
        return 0.0
    if mu == 0:
        return math.log(m / (m - nu)) / math.log(m / (m - u))
    B = m - u - mu
    return root({1: mu, B: 1 + (u - nu) / B}, m, 1)


rows = [('E8 P112 control', 9, 4 / 7), ('E8 P120', 9, 3 / 5), ('LZD_9  (h=9)', 9, 0.75)]
for n in list(range(11, 32, 2)) + [41, 51, 65]:
    rows.append((f'LZD_{n}', n, 6 / (n - 1)))
# FW_{n,k} exact values from scale.out (h, gamma)
rows += [('FW_{24,6}', 23, 48852 / 134596), ('FW_{32,6}', 31, 159216 / 906192), ('FW_{48,6}', 47, 831336 / 12271512),
         ('FW_{24,8}', 23, 250746 / 735471), ('FW_{32,8}', 31, 1130136 / 10518300), ('FW_{48,8}', 47, 9198276 / 377348994),
         ('gamma=0, h=9 (free)', 9, 0.0), ('gamma=0, h=13', 13, 0.0), ('gamma=0, h=31', 31, 0.0)]
ref = {}
print(f"{'system':16s} {'h':>3s} {'gamma':>6s} | lp(data floor) | selected-word s=9 | peel mu=0 s=9 (a) | peel mu=1 s=9 (a)")
for nm, h, g in rows:
    kk, s = lp(h, 1, g)
    vals = [kk, sel(h, g), peel(h, g), peel(h, g, mu=1)]
    if nm.startswith('E8 P112'):
        ref = vals
    if nm.startswith('gamma=0, h=9'):
        free9 = [v / ref[i] for i, v in enumerate(vals)]
    cells = [f"{v:.5f} (x{v/ref[i]:.2f} vs E8)" for i, v in enumerate(vals)]
    print(f"{nm:16s} {h:3d} {g:6.3f} | " + " | ".join(cells) + (f" s{s}" if s else ''))
print("lp/selected columns are kappa ceilings; peel columns are the tilt root a (kappa ~ a/(1+a) unfused).")
print("E8 control absolute: lp x%.2f, selected x%.2f over 0.0017" % (ref[0] / 0.0017, ref[1] / 0.0017))
assert abs(free9[0] - 5.02) < 0.01 and abs(free9[1] - 2.80) < 0.01, free9   # the u = 1 decoder-lane ceiling
print('PASS_EXPXA_PRICE')
