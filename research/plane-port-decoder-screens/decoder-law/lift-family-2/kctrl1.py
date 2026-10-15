"""Exact (Fractions) control of the |S| = 1 Y-line classification used by s1relax.py / s1gs.py / hz.py.
Builds kappa-family members A = {s}, B = supp kap' subset C(s), with structured level functions L, computes ALL Y-deficient
lines by exact kernels over the 510 classes, and checks:
 (1) every deficient line is a pencil line (u ~ r_s) or generic with Def subset {H0, h_b, h_b*} (b with u0 ~ p_b);
 (2) pencil lines: Def = {h : C-part(h) subset of the level block}  (dimD);
 (3) generic line with h_b exists  <=>  L const on M_b^- and L_b differs;
 (4) H0 lies on a generic line  <=>  L affine non-constant on N;  kernels never exceed dim 1;
 (5) X side: exactly e9 (Def q_s^perp) and g_s (Def <q_s>);
 (6) end-to-end: maxdef.mincost (exact-rho gauge floor over all dual-basis pairs of these lines) >= 72.   Usage: kctrl1.py seed ntrials"""
import sys, random, itertools
from fractions import Fraction as Fr
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from kfam import RT, dot, kernel, normline, CLASS
from s1lib import *
rnd = random.Random(int(sys.argv[1]))
def rq(): return Fr(rnd.choice([-1, 1]) * rnd.randint(1, 30), rnd.randint(1, 30))
RS = RT[s]
sgn = {b: (1 if dot(RT[b], RS) > 0 else -1) for b in N}
pF = {b: tuple(sgn[b] * x - y / 2 for x, y in zip(RT[b], RS)) for b in N}
def build(Lev, Ls_half):
    """Lev: dict b -> level (b in N); Ls_half: level of s (None = infinity, alp'_s = 0)."""
    al_s, ka_s = rq(), rq()
    X = []; Y = []
    for a in range(NL):
        if a == s:
            if Ls_half is None:
                als_p = Fr(0); kas_p = 1 / ka_s
            else:
                # level(s) = kap'_s / (2 alp'_s) = Ls_half and 2 al alp' + ka kap' = 1
                als_p = 1 / (2 * al_s + 2 * ka_s * Ls_half); kas_p = 2 * Ls_half * als_p
            X.append([al_s * c for c in RT[a]] + [ka_s]); Y.append([als_p * c for c in RT[a]] + [kas_p])
        else:
            X.append(list(RT[a]) + [Fr(0)])
            kp = Fr(0) if a in E else Lev[a] * Fr(1, 2) * sgn[a]       # L_b = kap'_b/(alp'_b r_b.r_s), alp' = 1/2
            Y.append([c / 2 for c in RT[a]] + [kp])
    for a in range(NL):
        for b in range(NL):
            if a == b or IP[a][b] == 0:
                assert dot(X[a], Y[b]) == (1 if a == b else 0)
    return X, Y
def deflines(M):
    out = {}
    for h, O in CLASS.items():
        Kk = kernel([M[a] for a in O])
        assert len(Kk) <= 1, ('kernel dim', h, len(Kk))
        if Kk: out.setdefault(normline(Kk[0]), []).append(h)
    return out
def ispencil(g):
    u = g[:8]; return all(u[i] * RS[j] == u[j] * RS[i] for i in range(8) for j in range(8))
def check(Lev, Ls_half, tag):
    X, Y = build(Lev, Ls_half)
    dX = deflines(X)
    e9 = tuple([Fr(0)] * 8 + [Fr(1)])
    ok5 = len(dX) == 2 and sorted(dX[e9]) == sorted(h for h in range(1, 511) if not dotq(h, LAB[s])) and \
        any(v == [LAB[s]] for v in dX.values())
    dY = deflines(Y)
    lev = dict(Lev); lev[s] = Ls_half if Ls_half is not None else 'inf'
    blocks = {}
    for c in C: blocks.setdefault(lev[c], []).append(c)
    predP = sorted(sorted(h for h in CL if CPART[h] <= set(bl)) for bl in blocks.values())
    predP = [x for x in predP if x]
    gotP = sorted(sorted(v) for g, v in dY.items() if ispencil(g))
    gen = [sorted(v) for g, v in dY.items() if not ispencil(g)]
    allowed = all(set(v) <= {H0, HB[b], HB[star[b]]} for v in gen for b in [next((b for b in N if HB[b] in v), N[0])])
    constM = {b for b in N if len({Lev[c] for c in Mm[b]}) == 1 and Lev[b] != Lev[Mm[b][0]]}
    gotB = {b for b in N if any(HB[b] in v for v in gen)}
    # affine test: L = w.p + c0 on N (solve least squares exactly via rank)
    rows = [list(pF[b]) + [Fr(1), Lev[b]] for b in N]
    import numpy as np
    A = np.array([[float(x) for x in r] for r in rows]); aff = np.linalg.matrix_rank(A, 1e-9) == np.linalg.matrix_rank(A[:, :-1], 1e-9)
    nonconst = len({Lev[b] for b in N}) > 1
    gotH0 = any(H0 in v for v in gen)
    import numpy as np
    from maxdef import mincost
    tov = lambda g: (lambda v: v / np.linalg.norm(v))(np.array([float(x) for x in g]))
    mc = mincost([(tov(g), v) for g, v in dX.items()], [(tov(g), v) for g, v in dY.items()])[0]
    ok = mc >= 72 and ok5 and predP == gotP and allowed and constM == gotB and gotH0 == (aff and nonconst)
    print(tag, 'Xok', ok5, 'pencil ok', predP == gotP, 'generic allowed', allowed, 'h_b lines', len(gotB), '== pred', constM == gotB,
          'H0 gen', gotH0, 'aff', aff and nonconst, 'Ylines', len(dY), 'mincost', mc, 'OK' if ok else 'FAIL', flush=True)
    return ok
CLQ = None
def clique7():
    best = []
    def ext(Q, cand):
        nonlocal best
        if len(Q) > len(best): best = Q
        for c in cand:
            if not Q or c > Q[-1]: ext(Q + [c], [d for d in cand if d > c and PP[c, d] == 1])
    ext([], N); return best
Q7 = clique7()
allok = True; n = 0
for tr in range(int(sys.argv[2])):
    kind = tr % 6
    v = rq()
    if kind == 0:                     # generic levels
        Lev = {b: rq() for b in N}; Ls = rq()
    elif kind == 1:                   # Q-config: constant block containing all M_b^-, b in Q at other levels
        Q = rnd.sample(Q7, rnd.randint(1, 7)); Lev = {b: v for b in N}
        common = rnd.random() < 0.5; w = rq()
        for b in Q: Lev[b] = w if common else rq()
        Ls = rnd.choice([v, rq(), None])
    elif kind == 2:                   # antipodal blocks
        b = rnd.choice(N); Lev = {c: rq() for c in N}
        v2 = rq()
        for c in Mm[b]: Lev[c] = v
        for c in Mp[b]: Lev[c] = v2
        Ls = rnd.choice([v, v2, rq(), None])
    elif kind == 3:                   # affine along p_b
        b = rnd.choice(N); c0 = rq(); B_ = rq()
        Lev = {c: c0 + B_ * PP[b, c] for c in N}; Ls = rnd.choice([c0 + B_, rq(), None])
    elif kind == 4:                   # affine along a random root of E7 (w = r_e, e in E)
        e = rnd.choice(E); c0 = rq(); B_ = rq()
        Lev = {c: c0 + B_ * dot(pF[c], RT[e]) for c in N}; Ls = rnd.choice([c0, rq()])
    else:                             # big block minus one line (engine extremal) and s inside
        b = rnd.choice(N); Lev = {c: v for c in N}; Lev[b] = rq(); Ls = v
    if kind != 5 and kind != 1 and len({Lev[c] for c in N}) == 1: continue
    ok = check(Lev, Ls, f'trial {tr} kind {kind}'); allok &= ok; n += 1
print('PASS' if allok else 'FAIL', n)
