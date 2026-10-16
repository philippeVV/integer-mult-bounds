"""a2credit-style exact label enumeration, generic. Xl/Yl: name -> Def subspace (set of F2^9 ints incl 0).
forced: set of (x, y) pairs whose pairing is nonzero (they must share a position if both used);
zero: set of (x, y) pairs whose pairing is 0 (cannot share a position). All other pairs free (relaxation).
indep_X / indep_Y: optional predicates on a tuple of names (linear independence); default: no constraint (relaxation)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys
pass   # expSM modules (e8, kfam, maxdef, cliques_lib, orbits, catalog, subsys) are copied into this directory
from maxdef import rho_def
def max_credit(Xl, Yl, forced, zero=frozenset(), indep_X=None, indep_Y=None):
    xs = list(Xl); ys = list(Yl); best = [0, None]
    c0 = {('x', x): 9 - rho_def(Xl[x], {0}) for x in xs}; c0.update({('y', y): 9 - rho_def({0}, Yl[y]) for y in ys})
    cp = {(x, y): 9 - rho_def(Xl[x], Yl[y]) for x in xs for y in ys}
    def rec_x(i, XS):
        if i == len(xs):
            if indep_X and not indep_X(tuple(XS)): return
            rec_y(0, XS, []); return
        rec_x(i + 1, XS); rec_x(i + 1, XS + [xs[i]])
    def rec_y(j, XS, YS):
        if j == len(ys):
            if indep_Y and not indep_Y(tuple(YS)): return
            match(XS, YS); return
        rec_y(j + 1, XS, YS); rec_y(j + 1, XS, YS + [ys[j]])
    def match(XS, YS):
        def rec(i, used, M, cr):
            if i == len(XS):
                Mx = dict(M)
                for x in XS:
                    for y in YS:
                        if (x, y) in forced and Mx.get(x) != y: return
                pos = len(XS) + len(YS) - len(M)
                if pos > 9: return
                tot = cr + sum(c0[('y', y)] for y in YS if y not in used)
                if tot > best[0]: best[0] = tot; best[1] = (list(M), XS, YS)
                return
            x = XS[i]
            rec(i + 1, used, M, cr + c0[('x', x)])
            for y in YS:
                if y not in used and (x, y) not in zero: rec(i + 1, used | {y}, M + [(x, y)], cr + cp[(x, y)])
        rec(0, frozenset(), [], 0)
    rec_x(0, [])
    return best
