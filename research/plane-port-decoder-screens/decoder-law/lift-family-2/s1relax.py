"""CP-SAT upper bound on the Y-side total deficiency T_Y for |S| = 1 (any B subset C(s), any levels), e9 not used.
Labels: 1, 2 = the (at most two) pencil blocks in the decoder; 3, 4 = other blocks holding a constant M_b^- (<= 2 such,
since non-antipodal M's intersect: facts.py); 0 = everything else (no constancy claimed).
Valid relaxation: every true configuration maps to a feasible assignment with >= objective."""
import sys
from ortools.sat.python import cp_model
from s1lib import *
def build(opts):
    m = cp_model.CpModel()
    K = [0, 1, 2, 3, 4]
    lab = {(c, k): m.NewBoolVar('') for c in C for k in K}
    for c in C: m.AddExactlyOne(lab[c, k] for k in K)
    d = {}; used = {}
    for k in (1, 2):
        if not opts.get('norank'): m.Add(sum(lab[c, k] for c in C) <= len(C) - 1)        # Y rank 9
        inD = {}
        for h, cp in CPART.items():
            v = m.NewBoolVar(''); inD[h] = v
            for c in cp: m.AddImplication(v, lab[c, k])
        e = [m.NewBoolVar('') for j in range(1, 10)]
        for j in range(1, 10):
            m.Add(sum(inD.values()) >= (2 ** j - 1)).OnlyEnforceIf(e[j - 1])
            if j > 1: m.AddImplication(e[j - 1], e[j - 2])
        d[k] = sum(e); used[k] = e[0]
    gb = {}
    for b in N:
        g = m.NewBoolVar(''); gb[b] = g; zs = []
        for k in (1, 2, 3, 4):
            z = m.NewBoolVar(''); zs.append(z)
            for c in Mm[b]: m.AddImplication(z, lab[c, k])
            if not opts.get('nob'): m.AddImplication(z, lab[b, k].Not())
        m.Add(sum(zs) == g)
    for b in N:
        if b < star[b] and not opts.get('nodim'): m.Add(gb[b] + gb[star[b]] + used[1] + used[2] <= 3)
    m.Add(used[1] + used[2] + sum(gb.values()) <= 9)
    for f in opts.get('fix', []): f(m, lab, gb, d, used)
    m.Maximize(d[1] + d[2] + sum(gb.values()))
    return m, lab, gb, d, used
def solve(opts={}, tl=600):
    m, lab, gb, d, used = build(opts)
    sv = cp_model.CpSolver(); sv.parameters.num_workers = 1; sv.parameters.max_time_in_seconds = tl
    r = sv.Solve(m)
    st = sv.StatusName(r)
    blocks = {k: [c for c in C if sv.Value(lab[c, k])] for k in (1, 2, 3, 4)} if r in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    return st, sv.ObjectiveValue(), sv.BestObjectiveBound(), blocks, ([b for b in N if sv.Value(gb[b])] if blocks else None), \
        ((sv.Value(d[1]), sv.Value(d[2])) if blocks else None)
if __name__ == '__main__':
    o = {a: True for a in sys.argv[1:]}
    st, val, bd, blocks, Q, dd = solve(o); print('opts', o)
    print(st, 'max T_Y (no h0-generic)', val, 'bound', bd, 'pencil dims', dd, 'Q', Q)
    print({k: len(v) for k, v in blocks.items()}, 's in', [k for k, v in blocks.items() if s in v])
    print('recheck dims', dimD(blocks[1]), dimD(blocks[2]))
