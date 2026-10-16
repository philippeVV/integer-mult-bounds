"""Per R3 orbit / structured dependent non-clique support A:
 (1) row-wise exact builder: random concurrent c -> feasible, Y rows exactly kappa-form (8-part ~ r_a), freedom only on F;
     one non-concurrent perturbation of a single c_s -> infeasible.
 (2) exact kappa-family members (A, random B subset C(A), alpha-zero patterns): exact deficient lines == Lemma 2 prediction,
     label enumeration (credit.py, all non-forced pairings free) <= 9, exact-rho gauge floor mincost >= 72.
Usage: members.py seed n_per_support [ctrl]  (ctrl: drop the forced e9.h'_c, g_f.e9' pairings; label credit must exceed 9)"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, json, random, itertools
from r3lib import *
from kfam import member
from credit import max_credit
from maxdef import span2
from cliques_lib import common
from orbits import clean
rnd = random.Random(int(sys.argv[1])); NPER = int(sys.argv[2]); CTRL = sys.argv[3:] == ['ctrl']
def rq(): return Fr(rnd.choice([-1, 1]) * rnd.randint(1, 9), rnd.randint(1, 9))
def span_lines(S):
    u, s, vt = np.linalg.svd(RF[list(S)]); Q = vt[:(s > 1e-9).sum()]
    res = (RF * RF).sum(1) - ((RF @ Q.T) ** 2).sum(1); return [a for a in range(NL) if res[a] < 1e-9]
sup = [('orb4 ' + str(o['rep']), o['rep']) for o in json.load(open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'orb4.json')))]
A3 = span_lines(sup[2][1]); sup.append(('A3 (6)', A3))
D4 = None
for q in itertools.combinations(range(40), 4):
    L = span_lines(q)
    if len(L) == 12 and rank_S(L) == 4: D4 = L; break
sup += [('D4 (12)', D4), ('D4 9 lines', D4[:9]), ('D4 circuit 5', None)]
# D4 circuit: 4 mutually orth lines + 1 line non-orth to all of them inside D4
for q in itertools.combinations(D4, 4):
    if all(IP[a][b] == 0 for a, b in itertools.combinations(q, 2)):
        e = next((b for b in D4 if all(IP[b][a] != 0 for a in q)), None)
        if e is not None: sup[-1] = ('D4 circuit 5 (4 orth + 1)', list(q) + [e]); break
while True:  # a random dependent non-clique 20-set
    S = rnd.sample(range(NL), 20)
    if max(sum(1 for s in S if IP[a][s] == 0 and s != a) for a in range(NL)) <= 26: break
sup.append(('random 20', S))
def predict(A, P, side):
    L = {}; e9 = tuple([Fr(0)] * 8 + [Fr(1)])
    D = [h for h in range(1, 511) if all(dotq(LAB[a], h) == 0 for a in A)]
    if D: L[e9] = sorted(D)
    for f in A:
        x, k, y, kp = P[f]; al, ka = (x, k) if side == 'X' else (y, kp)
        if all(IP[f][a] != 0 for a in A if a != f) or al == 0:
            g = normline([ka * c for c in RT[f]] + [-2 * al]); L.setdefault(g, []).append(LAB[f]); L[g] = sorted(L[g])
    return L
allok = True; worst = 99
for name, A in sup:
    A = sorted(A); C = common(A); T = [t for t in range(NL) if t not in A]
    F = [t for t in T if all(IP[t][s] != 0 for s in A)]
    # (1) builder
    c0 = [rq() for _ in range(8)]; mu = {s: rq() for s in A}
    c = {s: [x + mu[s] * y for x, y in zip(c0, RT[s])] for s in A}; xi = {s: rq() for s in A}
    out = build(A, c, xi, rnd)
    okb = out is not None
    if okb:
        X, Y, fr = out
        okb = check(X, Y) and {a for a in range(NL) if fr[a]} == set(C) and all(fr[a] <= 1 for a in range(NL))
        # kappa form after the c0 shift: y_a 8-part ~ r_a
        okb &= all(len(kernel([Y[a][:8], list(RT[a])], 8)) == 7 for a in range(NL))
    s0 = A[0]; c2 = dict(c); c2[s0] = [x + rq() for x in c[s0]]
    infeas = build(A, c2, xi, rnd) is None
    rkY = np.linalg.matrix_rank(np.array([[float(x) for x in y] for y in Y])) if out else None
    if not C:
        print(f'{name:34s} |A|={len(A):2d} C(A) empty: builder-kappa-form={okb} rank Y={rkY} (=8: no rank-9 lift on this support) nonconc-infeasible={infeas}', flush=True)
        allok &= okb and infeas and rkY == 8; continue
    # (2) kappa members
    res = []
    for tr in range(NPER):
        mode = tr % 3
        B = [b for b in C if rnd.random() < (0.5 if mode == 0 else 0.15)] if mode < 2 else rnd.sample(C, min(len(C), rnd.randint(1, 3)))
        if not B: B = [rnd.choice(C)]
        AB = [a for a in A if a in B]
        al = {a: Fr(0) for a in AB if rnd.random() < 0.4}; alp = {a: Fr(0) for a in AB if a not in al and rnd.random() < 0.4}
        X, Y, P = member(A, B, rnd, al, alp); assert check(X, Y)
        mc, dX, dY, kd = gauge_floor(X, Y)
        cleanAB = not clean(sum(1 << a for a in A)) and not clean(sum(1 << b for b in B))
        okp = ({l: sorted(h) for l, h in dX.items()} == predict(A, P, 'X') and {l: sorted(h) for l, h in dY.items()} == predict(B, P, 'Y'))
        # label enumeration on the exact lines; forced = Theorem-1 pairings (every nonzero pairing with e9 or e9', and g_f.h'_f); g_f.h'_c (f != c) free
        e9 = tuple([Fr(0)] * 8 + [Fr(1)])
        Xl = {('X', l): span2(h) for l, h in dX.items()}; Yl = {('Y', l): span2(h) for l, h in dY.items()}
        forced = {(x, y) for x in Xl for y in Yl if dot(x[1], y[1]) != 0 and ((x[1] == e9 or y[1] == e9) and not CTRL or x[1] == e9 and y[1] == e9 or
                  any(normline([P[f][1] * v for v in RT[f]] + [-2 * P[f][0]]) == x[1] and normline([P[f][3] * v for v in RT[f]] + [-2 * P[f][2]]) == y[1] for f in AB))}
        zero = {(x, y) for x in Xl for y in Yl if dot(x[1], y[1]) == 0}
        cr = max_credit(Xl, Yl, forced, zero)[0]
        ok = kd == (1, 1) and cleanAB and okp and (CTRL or cr <= 9) and mc is not None and mc >= 72
        res.append((len(B), len(al) + len(alp), cr, mc, ok)); allok &= ok; worst = min(worst, mc if mc is not None else -1)
    print(f'{name:34s} |A|={len(A):2d} |C(A)|={len(C):2d} |F|={len(F):2d} builder-kappa-form={okb} nonconc-infeasible={infeas} | members: credit max {max(r[2] for r in res)}, mincost {sorted({r[3] for r in res})}, all OK {all(r[4] for r in res)}', flush=True)
    allok &= okb and infeas and (not CTRL or max(r[2] for r in res) > 9)
print('PASS' if allok else 'FAIL', 'min mincost', worst)
print('PASS_EXPR3_MEMBERS' if allok else 'FAIL_EXPR3_MEMBERS')
