"""|S| = 2 lifts with a NON-Clean B (Theorem 1 does not apply): B contains F up to <= 1 line. Relaxed credit bound with the
extra (non-standard) Y-lines included, each with an UPPER bound on its Def (credit is monotone in Def):
  Dmax(u) = span{h : O_h minus B lies in u^perp}, for u' in the A2 plane P = span(r_s1, r_s2) (special directions r_s1, r_s2,
  r_t (third A2 line) and generic), all extra lines have gamma' != 0 (every non-Clean class contains all of F, rank R[F] = 8; checked)
  -> forced nonzero pairing with e9; pairing with g_f free. All lines of the form (u' in P, gamma') and e9', g_s, e9 lie in the
  3-dim space W3 = P x K: at most 3 Y-lines and 3 X-lines from W3.  Zero patterns (al/alp = 0 on A cap B) chosen by CP-SAT.
Usage: relax2.py  -> enumerates all non-Clean B (F minus <= 1 line, any S cap B)."""
import json, itertools, numpy as np
from ortools.sat.python import cp_model
from e8 import *
from maxdef import rho_def, span2
from orbits import clean, CL, R
orb = json.load(open('orbits.json')); S = orb['2'][0]['rep']; C = orb['2'][0]['C']; F = [c for c in C if c not in S]
s1, s2 = S; t = next(c for c in F if np.linalg.matrix_rank(R[[s1, s2, c]]) == 2)
Fm = sum(1 << b for b in F); Cm = sum(1 << c for c in C)
BAD = [h for h, _ in clean(Cm)]
assert np.linalg.matrix_rank(R[F]) == 8 and all(CL[h] & Fm == Fm for h in BAD), 'gamma=0 extra lines impossible'
assert not clean(sum(1 << a for a in S)), 'X side (A = S) must be Clean'
print('Clean(S) holds (X side standard);', 'bad classes of C(S):', BAD, '; each contains F, rank R[F] = 8  => every extra line has gamma != 0')
def subspace(hs):
    return span2(sorted(set(hs)))
def dset(u, Bm):
    return [h for h, O in CL.items() if all(abs(R[a] @ u) < 1e-9 for a in bits(O & ~Bm))]
def dmax(u, Bm): return subspace(dset(u, Bm))
def h_nonzero(D): return len(D) > 1
def run(B, ZXa, ZYa):
    Bm = sum(1 << b for b in B); A = S
    UA = [LAB[a] for a in A]; UB = [LAB[b] for b in B]
    perp = lambda labs: {h for h in range(512) if all(dotq(l, h) == 0 for l in labs)}
    X = [('e9', perp(UA), True)] + [(f, dmax(R[f], sum(1 << a for a in A)) if False else {0, LAB[f]}, True) for f in A]
    Y = [("e9'", perp(UB), True)]
    for c in B:
        if all(IP[c][b] != 0 for b in B if b != c) or c in S:
            Y.append((('h', c), dmax(R[c], Bm), c in S or c == t))          # gamma' != 0 version (merged allowed)
            if c in S: Y.append((('z', c), {0, LAB[c]}, True))                # gamma' = 0 version (alp_c = 0): Def exactly <q_c>
    P2 = np.array([R[s1], R[s2]]); n = lambda r: r - (P2.T @ np.linalg.solve(P2 @ P2.T, P2 @ r))  # unused helper
    def perpP(r):   # direction in P orthogonal to r (r in P)
        a, b = R[s1], R[s2]; M = np.array([[a @ a, a @ b], [a @ b, b @ b]]); c = np.linalg.solve(M, [0, 1]) if abs(r @ a) > 1e-9 else None
        v = a * (r @ b) - b * (r @ a); return v
    dirs = [('E1', R[s1], s1), ('E2', R[s2], s2), ('Et', R[t], None), ('N1', perpP(R[s1]), None), ('N2', perpP(R[s2]), None),
            ('Nt', perpP(R[t]), None), ('EPg', R[s1] + 0.37 * R[s2], None)]
    for name, u, si in dirs:
        st = dset(u, Bm)
        if si is not None: st = [h for h in st if not (CL[h] >> si) & 1]     # lines of Pi_si other than h'_si
        for cp in range(2):
            Y.append((f'E:{name}:{cp}:full', subspace(st), True))
            if 499 in st: Y.append((f'E:{name}:{cp}:no499', subspace([h for h in st if h != 499]), True))
    for c in list(B):
        if all(IP[c][b] != 0 for b in B if b != c) or c in S:
            st = dset(R[c], Bm)
            if 499 in st and c in S: Y.append((('h', c, 'no499'), subspace([h for h in st if h != 499]), True))
    return X, Y
def solve(B):
    X, Y = run(B, None, None); AB = [a for a in S if a in B]
    m = cp_model.CpModel(); zx = {a: m.NewBoolVar('zx%d' % a) for a in AB}; zy = {a: m.NewBoolVar('zy%d' % a) for a in AB}
    ux = [m.NewBoolVar('') for _ in X]; uy = [m.NewBoolVar('') for _ in Y]; obj = []
    cx = [9 - rho_def(D, {0}) for _, D, _ in X]; cy = [9 - rho_def({0}, D) for _, D, _ in Y]
    Mv = {}
    for i, (x, Dx, _) in enumerate(X):
        for j, (y, Dy, _) in enumerate(Y):
            v = m.NewBoolVar(''); Mv[i, j] = v; m.AddImplication(v, ux[i]); m.AddImplication(v, uy[j])
            both = [ux[i], uy[j]]
            if x == 'e9' and y == "e9'": m.AddBoolOr([l.Not() for l in both] + [v])
            elif x == 'e9' and isinstance(y, tuple) and y[0] == 'z': pass
            elif x == 'e9' and isinstance(y, tuple):
                c = y[1]
                if c in zy: m.AddBoolOr([l.Not() for l in both] + [v, zy[c]]); m.AddImplication(zy[c], v.Not())
                else: m.AddBoolOr([l.Not() for l in both] + [v])
            elif x == 'e9': m.AddBoolOr([l.Not() for l in both] + [v])            # extra lines: gamma' != 0
            elif y == "e9'":
                if x in zx: m.AddBoolOr([l.Not() for l in both] + [v, zx[x]]); m.AddImplication(zx[x], v.Not())
                else: m.AddBoolOr([l.Not() for l in both] + [v])
            elif isinstance(y, tuple) and y[0] == 'z': pass
            elif isinstance(y, tuple) and y[1] == x: m.AddBoolOr([l.Not() for l in both] + [v])
            elif isinstance(y, tuple):
                c = y[1]; lits = [zx[x]] if x in zx else []
                if c in zy: lits.append(zy[c])
                for z in lits: m.AddBoolOr([l.Not() for l in both] + [v, z.Not()])
            obj.append(-(cx[i] + cy[j] - (9 - rho_def(Dx, Dy))) * v)
    for j, (y, D, _) in enumerate(Y):
        if isinstance(y, tuple) and y[0] == 'z':
            if y[1] in zy: m.AddImplication(uy[j], zy[y[1]])
            else: m.Add(uy[j] == 0)
            for i, (x, _, _) in enumerate(X):          # (r_c, 0): pairs 0 with e9, nonzero (kappa_f r_f.r_c) with every g_f
                if x == 'e9': m.Add(Mv[i, j] == 0)
                else: m.AddBoolOr([ux[i].Not(), uy[j].Not(), Mv[i, j]])
        if isinstance(y, tuple) and y[0] == 'h' and y[1] in zy: m.AddImplication(uy[j], zy[y[1]].Not())
    isno = lambda y: (isinstance(y, tuple) and y[-1] == 'no499') or (isinstance(y, str) and y.endswith('no499'))
    full499 = [j for j, (y, D, _) in enumerate(Y) if 499 in D and not isno(y)]
    m.Add(sum(uy[j] for j in full499) <= 2)
    groups = {}
    for j, (y, D, _) in enumerate(Y):
        key = y[:2] if isinstance(y, tuple) and len(y) == 3 else (y.rsplit(':', 1)[0] if isinstance(y, str) and y.startswith('E:') else y)
        groups.setdefault(str(key), []).append(j)
    for js in groups.values(): m.Add(sum(uy[j] for j in js) <= 1)
    for i in range(len(X)): m.Add(sum(Mv[i, j] for j in range(len(Y))) <= 1)
    for j in range(len(Y)): m.Add(sum(Mv[i, j] for i in range(len(X))) <= 1)
    m.Add(sum(ux) + sum(uy) - sum(Mv.values()) <= 9)
    m.Add(sum(uy[j] for j, (_, _, w) in enumerate(Y) if w) <= 3)
    m.Maximize(sum(a * b for a, b in zip(cx, ux)) + sum(a * b for a, b in zip(cy, uy)) + sum(obj))
    s = cp_model.CpSolver(); s.parameters.num_workers = 1; assert s.Solve(m) == cp_model.OPTIMAL
    conf = ([X[i][0] for i in range(len(X)) if s.Value(ux[i])], [Y[j][0] for j in range(len(Y)) if s.Value(uy[j])])
    return int(round(s.ObjectiveValue())), conf, [len(D).bit_length() - 1 for _, D, _ in Y]
if __name__ == '__main__':
    mx = 0; n = 0
    for drop in [None] + F:
        for SB in itertools.chain.from_iterable(itertools.combinations(S, r) for r in range(3)):
            B = [b for b in F if b != drop] + list(SB)
            if not clean(sum(1 << b for b in B)): continue
            v, conf, dims = solve(B); n += 1; mx = max(mx, v)
            if v > 9 or n <= 3: print('B = F' + ('' if drop is None else f' - {drop}'), '+', SB, 'credit', v, conf, 'Y Def dims', dims, flush=True)
    print('non-Clean B tested', n, 'max credit', mx, 'PASS' if mx <= 9 else 'CREDIT>9 (relaxation inconclusive)')
