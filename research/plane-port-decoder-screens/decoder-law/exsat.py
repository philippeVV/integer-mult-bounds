"""Exhaustive Boolean-support relaxation over a CLASS of support pairs on P120 (or any L):
choose any multiset of candidate support pairs (R_c, A_c) with total pairing rank <= BUDGET such that the
Boolean relaxation of an exact masked decoder is satisfiable.  UNSAT => no exact decoder (over any field)
whose terms all have supports inside the class with cost <= BUDGET.
Relaxation (necessary): x_ac = [v_c(a) != 0] (a in R_c), y_bc = [w_c(b) != 0] (b in A_c); every diagonal needs a
nonzero product; a required zero cannot have exactly one nonzero product."""
import sys, time
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from e8 import *
def run(L, cands, budget, verbose=True):
    """cands: list of (R, A, cost, mult). Returns (sat, chosen list)."""
    t0 = time.time(); pool = IDPool(); s = Cadical153(); ncl = 0
    def add(c):
        nonlocal ncl; s.add_clause(c); ncl += 1
    slots = []   # (cand index, copy)
    for ci, (R, A, c, m) in enumerate(cands):
        for j in range(m): slots.append((ci, j))
    sel = {sl: pool.id(('s',) + sl) for sl in slots}
    X = {}; Y = {}
    for sl in slots:
        R, A, c, m = cands[sl[0]]
        for a in bits(R & L): X[a, sl] = pool.id(('x', a) + sl); add([-X[a, sl], sel[sl]])
        for b in bits(A & L): Y[b, sl] = pool.id(('y', b) + sl); add([-Y[b, sl], sel[sl]])
    # symmetry breaking between copies: copy j used => copy j-1 used
    for (ci, j) in slots:
        if j > 0: add([-sel[(ci, j)], sel[(ci, j - 1)]])
    Ls = list(bits(L))
    for a in Ls:
        for b in Ls:
            if a != b and IP[a][b] != 0: continue
            ps = []
            for sl in slots:
                if (a, sl) in X and (b, sl) in Y:
                    p = pool.id(('p', a, b) + sl); x, y = X[a, sl], Y[b, sl]
                    add([-p, x]); add([-p, y]); add([p, -x, -y]); ps.append(p)
            if a == b:
                if not ps: return False, None
                add(ps); continue
            if len(ps) < 2:
                if ps: add([-ps[0]])
                continue
            # not exactly one: p_i -> (OR_{j<i} p_j) or (OR_{j>i} p_j)
            m = len(ps); pre = [None] * m; suf = [None] * m
            for i in range(1, m):
                pre[i] = pool.id(('pre', a, b, i)); add([-pre[i], ps[i - 1]] + ([pre[i - 1]] if i > 1 else []))
            for i in range(m - 2, -1, -1):
                suf[i] = pool.id(('suf', a, b, i)); add([-suf[i], ps[i + 1]] + ([suf[i + 1]] if i < m - 2 else []))
            for i in range(m):
                add([-ps[i]] + ([pre[i]] if pre[i] else []) + ([suf[i]] if suf[i] is not None else []))
    lits = []
    for sl in slots: lits += [sel[sl]] * 0 + [pool.id(('w', k) + sl) for k in range(cands[sl[0]][2])]
    for sl in slots:
        for k in range(cands[sl[0]][2]): add([-sel[sl], pool.id(('w', k) + sl)])
    enc = CardEnc.atmost(lits=lits, bound=budget, top_id=pool.top, encoding=EncType.seqcounter)
    for c in enc.clauses: add(c)
    if verbose: print(f'  built: slots {len(slots)} clauses {ncl} vars {pool.top} ({time.time() - t0:.0f}s)', flush=True)
    r = s.solve()
    chosen = None
    if r:
        mdl = set(l for l in s.get_model() if l > 0)
        chosen = [sl for sl in slots if sel[sl] in mdl]
    s.delete()
    if verbose: print(f'  solve: {"SAT" if r else "UNSAT"} ({time.time() - t0:.0f}s)', flush=True)
    return r, chosen
HYP = {h: sum(1 << a for a in range(NL) if dotq(LAB[a], h) == 0) for h in range(1, 512)}
if __name__ == '__main__':
    mode, budget = sys.argv[1], int(sys.argv[2])
    L = FULL
    if mode == 'ctrl':   # the 8 incumbent frame supports only (cost 72): budget 72 SAT, 71 UNSAT
        F = frames()[0]; cands = [(supp(ROOT[f], FULL), supp(ROOT[f], FULL), 9, 1) for f in F]
    elif mode == 'symhyp':   # all symmetric hyperplane-section supports (+ full support), 1 copy each hyperplane
        cands = [(m, m, rho(m, m), int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3].isdigit() else 1) for h, m in HYP.items() if m and rho(m, m) <= 8]
        if "nofull" not in sys.argv: cands.append((FULL, FULL, 9, 8))
    elif mode == 'symhyp_rho7':
        cands = [(m, m, rho(m, m), int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3].isdigit() else 1) for h, m in HYP.items() if m and rho(m, m) <= 7]
        if "nofull" not in sys.argv: cands.append((FULL, FULL, 9, 8))
    elif mode == 'e64':   # only the 120 E7+A1 sections (rho 7, 64 lines)
        cands = [(m, m, 7, 1) for h, m in HYP.items() if bin(m).count('1') == 64 and rho(m, m) == 7]
    elif mode == 'e56':
        cands = [(m, m, 7, 1) for h, m in HYP.items() if bin(m).count('1') == 56 and rho(m, m) == 7]
    print(mode, 'budget', budget, 'candidates', len(cands), flush=True)
    r, ch = run(L, cands, budget)
    if ch:
        print('chosen', [(cands[c][2], bin(cands[c][0]).count('1')) for c, j in ch])
        inv = {m: h for h, m in HYP.items()}
        print('chosen_h', [(inv.get(cands[c][0]), inv.get(cands[c][1])) for c, j in ch], flush=True)
