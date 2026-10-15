"""Independent checker for LEMMA S (support SAT), k = 7 version: ../lemmaS/check_S.py with capacities up to 3
allowed (a group touching Q contributes e_g = min(c_g,|S_g&Q|) copies of its touch literal in (R)).
Imports nothing from the engine or model code.
Same as check_S7.py except that the static (R) slack bound (default 2) is read from the SLACK env var.

Setting: nine-term masked co-star decoder on P112 (E8 lines, PR #1320 order, frame F0 deleted).
Normals a_1..a_k (distinct), capacities c_g in {1,2,3}, sum 9.  Group g: D_g, T_g in K^{S_g x c_g},
S_g = live lines orthogonal to a_g or equal to a_g.  M = sum_g D_g T_g^T must satisfy M[r,r]=4 and
M[r,s]=0 for orthogonal live r,s (any field of char != 2).

Booleans  A[g,r] <-> D_g[r,:] != 0,   B[g,r] <-> T_g[r,:] != 0   (GL(c_g)-invariant).
 (D)  diagonal r: some g with A[g,r] and B[g,r]                         [M[r,r] = 4 != 0]
 (X)  orthogonal r,s and cap-1 group x with A[x,r], B[x,s]: some other group h with A[h,r], B[h,s]
      [else M[r,s] = u_x[r] v_x[s] != 0]
 (R)  orthogonal clique Q (static slack <= 2): sum over groups g touching Q (some A[g,.] and some B[g,.]
      on Q) of e_g = min(c_g, |S_g & Q|) is >= |Q|      [rank 4I_Q = |Q| <= sum rank D_g[Q] T_g[Q]^T]
Every decoder induces a satisfying assignment (auxiliaries defined one-directionally), so UNSAT on both
CaDiCaL 1.5.3 and Glucose 4 proves that no decoder with these normals/capacities exists.
usage: check_S7s.py cases.json [--frames]   (list of [normals, caps]); --frames: (R) only for the 1260 live
8-frames (Bessel: |S_a & Q| >= 4 = |Q| - 4, so e_g = c_g and slack = 9 - 8 = 1).
"""
import itertools, sys, json, hashlib
from pysat.formula import IDPool
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153, Glucose4

R8 = [y for y in itertools.product(range(-2, 3), repeat=8)
      if sum(v * v for v in y) == 8 and len({v % 2 for v in y}) == 1 and sum(y) % 4 == 0]
assert len(R8) == 240
def lead_pos(y):
    for v in y:
        if v: return y if v > 0 else tuple(-w for w in y)
LN = sorted({lead_pos(y) for y in R8},
            key=lambda v: (1, tuple(int(x < 0) for x in v[1:])) if all(abs(x) == 1 for x in v)
            else (0, tuple(i for i, x in enumerate(v) if x), int([x for x in v if x][1] < 0)))
assert hashlib.sha256(json.dumps([list(v) for v in LN], sort_keys=True, separators=(',', ':')).encode()).hexdigest() \
    == "1c71b1394110e20a71a73e50a4b7b23e3ffce25e8bd3c27e04a94b738fff798e"
dot = lambda a, b: sum(x * y for x, y in zip(LN[a], LN[b]))
DEL = [i for i in range(120) if sum(1 for x in LN[i] if x) == 2 and
       [j for j, x in enumerate(LN[i]) if x] in ([0, 1], [2, 3], [4, 5], [6, 7])]
assert DEL == [0, 1, 26, 27, 44, 45, 54, 55]
LV = [i for i in range(120) if i not in DEL]
PERP = {r: {s for s in LV if s != r and dot(r, s) == 0} for r in LV}
def costar(a): return {r for r in LV if r == a or dot(a, r) == 0}

# all orthogonal cliques via bitmasks
POS = {r: i for i, r in enumerate(LV)}
NB = [sum(1 << POS[s] for s in PERP[r]) for r in LV]
def all_cliques():
    out = []
    stack = [(1 << i, NB[i] & ~((2 << i) - 1)) for i in range(len(LV))]
    while stack:
        q, cand = stack.pop(); out.append(q)
        while cand:
            b = cand & -cand; i = b.bit_length() - 1; cand ^= b
            stack.append((q | b, cand & NB[i]))
    return [tuple(LV[i] for i in range(len(LV)) if q >> i & 1) for q in out]
CLQ = all_cliques()
assert len(CLQ) == 262928

def cnf(normals, caps):
    assert len(set(normals)) == len(normals) and sum(caps) == 9 and set(caps) <= {1, 2, 3}
    G = len(normals); Sg = [costar(a) for a in normals]
    vp = IDPool(); cl = []
    A = lambda g, r: vp.id(('A', g, r)); B = lambda g, r: vp.id(('B', g, r))
    def both(g, r, s):          # p -> A[g,r] & B[g,s]
        p = vp.id(('P', g, r, s))
        if p not in seen: seen.add(p); cl.append([-p, A(g, r)]); cl.append([-p, B(g, s)])
        return p
    seen = set()
    for r in LV:
        gs = [g for g in range(G) if r in Sg[g]]
        cl.append([both(g, r, r) for g in gs])
        for s in PERP[r]:
            hs = [g for g in gs if s in Sg[g]]
            for x in hs:
                if caps[x] == 1:
                    cl.append([-A(x, r), -B(x, s)] + [both(h, r, s) for h in hs if h != x])
    nR = 0
    for Q in CLQ:
        e = [min(caps[g], len(Sg[g].intersection(Q))) for g in range(G)]
        if sum(e) - len(Q) > SLACK or len(Q) < 2 or (FRAMES and len(Q) != 8): continue
        nR += 1
        lits = []
        for g in range(G):
            if e[g] == 0: continue
            t = vp.id(('T', g, Q)); rs = [r for r in Q if r in Sg[g]]
            cl.append([-t] + [A(g, r) for r in rs]); cl.append([-t] + [B(g, r) for r in rs])
            lits.append(t)
            for j in range(2, e[g] + 1):     # with t itself: e_g copies of the touch literal (rank <= e_g)
                tj = vp.id(('T%d' % j, g, Q)); cl.append([-tj, t]); lits.append(tj)
        cl += CardEnc.atleast(lits, len(Q), vpool=vp, encoding=EncType.seqcounter).clauses
    return cl, nR

FRAMES = "--frames" in sys.argv
import os; SLACK = int(os.environ.get("SLACK", "2"))
if __name__ == '__main__':
    cases = json.load(open(sys.argv[1])); bad = 0
    for normals, caps in cases:
        cl, nR = cnf(normals, caps)
        res = []
        for Sol in (Cadical153, Glucose4):
            s = Sol(bootstrap_with=cl); res.append(s.solve()); s.delete()
        bad += any(res)
        print(normals, caps, 'clauses', len(cl), 'rank-cliques', nR, 'cadical', 'SAT' if res[0] else 'UNSAT',
              'glucose', 'SAT' if res[1] else 'UNSAT', flush=True)
    print('LEMMA S:', 'all %d cases UNSAT on both solvers' % len(cases) if not bad else '%d cases NOT refuted' % bad)
