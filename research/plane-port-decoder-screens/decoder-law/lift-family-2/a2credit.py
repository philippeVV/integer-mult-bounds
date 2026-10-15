"""Credit bound for non-concurrent A2-triple lifts, all parameter values (README sec 2).
X lines (Lift Lemma, X_T rigid): e9 [Def span(q_S)^perp, dim 6], g_i [Def <q_i>], with g_3 in span(e9, g_1, g_2) => <= 3 of 4.
Y lines (a2scan*.py + proof): h'_i in K_i [Def <q_i>], unless merged h'_i = h'_j = L_ij [Def span(q_i, q_j)] (at most one merge).
Forced nonzero: g_i.h'_i (perfect pairing A_i x K_i, x_si in K_i, y_si in A_i, x_si.y_si = 1); e9.L_ij (gamma(L_ij) = 1).
All other pairings free (relaxation). Labels: q_i pairwise orthogonal, q_i.q_i = 1.
Usage: a2credit.py [ctrl]  (ctrl: drop the forced g_i.h'_i != 0; must exceed 9)"""
import sys, itertools
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from maxdef import rho_def, span2
from a2lib import S, LAB
ctrl = sys.argv[1:] == ['ctrl']
q = {i: LAB[x] for i, x in enumerate(S)}
De = {h for h in range(512) if all(bin(h & q[i]).count('1') % 2 == 0 for i in range(3))}
Xl = {'e9': De, 'g0': {0, q[0]}, 'g1': {0, q[1]}, 'g2': {0, q[2]}}
def ycases():
    yield 'none', {'h0': {0, q[0]}, 'h1': {0, q[1]}, 'h2': {0, q[2]}}, set()
    for i, j in itertools.combinations(range(3), 2):
        k = 3 - i - j
        yield f'M{i}{j}', {f'M{i}{j}': span2([q[i], q[j]]), f'h{k}': {0, q[k]}}, {(f'g{i}', f'M{i}{j}'), (f'g{j}', f'M{i}{j}'), ('e9', f'M{i}{j}')}
def forced(x, y, extra):
    if (x, y) in extra: return True
    if not ctrl and x[0] == 'g' and y == 'h' + x[1]: return True
    return False
best = (0, None)
for name, Yl, extra in ycases():
    xs = list(Xl); ys = list(Yl)
    for nx in range(len(xs) + 1):
        for XS in itertools.combinations(xs, nx):
            if len(XS) == 4: continue                       # rank 3
            for ny in range(len(ys) + 1):
                for YS in itertools.combinations(ys, ny):
                    # matchings between XS and YS (partial injective)
                    def rec(i, used, M):
                        global best
                        if i == len(XS):
                            Mx = dict(M); My = {v: k for k, v in M}
                            # unmatched cross pairs must pair 0: forbidden if forced nonzero
                            for x in XS:
                                for y in YS:
                                    if Mx.get(x) != y and forced(x, y, extra): return
                            pos = len(XS) + len(YS) - len(M)
                            if pos > 9: return
                            cr = 0
                            for x in XS:
                                if x in Mx: cr += 9 - rho_def(Xl[x], Yl[Mx[x]])
                                else: cr += 9 - rho_def(Xl[x], {0})
                            for y in YS:
                                if y not in My: cr += 9 - rho_def({0}, Yl[y])
                            if cr > best[0]: best = (cr, (name, M, XS, YS))
                            return
                        rec(i + 1, used, M)
                        for y in YS:
                            if y not in used: rec(i + 1, used | {y}, M + [(XS[i], y)])
                    rec(0, set(), [])
print('ctrl' if ctrl else 'true', 'max credit', best[0], best[1], 'PASS' if (best[0] <= 9) != ctrl else 'FAIL')
