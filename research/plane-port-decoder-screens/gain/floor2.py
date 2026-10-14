"""EXACT floor (any field):  C >= max over S (pairwise non-orthogonal planes) of |S| + |Q|, Q subset S with pairwise
rank-2 pairings (E_a x E_b of F2 rank 2).  Proof: N[S,S]=I_S; cost-1 terms touching SxS live on rank-1 pairs (off-diag),
so the cost>=2 terms carry I - X with X supported on rank-1 pairs, rank >= |Q|; #terms >= |S|.
CP-SAT finds S,Q; the witness is re-checked exactly and the optimality status is printed."""
import sys, json, itertools
from base import *
from ortools.sat.python import cp_model
def prk(Ea, Eb): return rank2([sum(dot(x, y) << j for j, y in enumerate(Eb)) for x in Ea])
def floor(r, t=60):
    E = planes(r)
    PR = [[prk(E[a], E[b]) for b in range(n)] for a in range(n)]
    m = cp_model.CpModel(); s = [m.NewBoolVar('') for _ in range(n)]; q = [m.NewBoolVar('') for _ in range(n)]
    for a, b in itertools.combinations(range(n), 2):
        if PR[a][b] == 0: m.AddBoolOr([s[a].Not(), s[b].Not()])
        if PR[a][b] != 2: m.AddBoolOr([q[a].Not(), q[b].Not()])
    for a in range(n): m.AddImplication(q[a], s[a])
    m.Maximize(sum(s) + sum(q)); sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = t; sv.parameters.num_workers = 4
    st = sv.Solve(m)
    S = [a for a in range(n) if sv.Value(s[a])]; Q = [a for a in range(n) if sv.Value(q[a])]
    assert all(PR[a][b] > 0 for a, b in itertools.combinations(S, 2)) and all(PR[a][b] == 2 for a, b in itertools.combinations(Q, 2)) and set(Q) <= set(S)
    return len(S), len(Q), sv.StatusName(st), int(sv.BestObjectiveBound())
for f in sys.argv[1:]:
    d = json.load(open(f)); runs = d if isinstance(d, list) else [d]
    for run in runs[:2]:
        r = run.get('r') or run.get('R')
        S, Q, st, ub = floor(r)
        print(f, '|S|', S, '|Q|', Q, 'floor C >=', S + Q, st, 'ub', ub, '(old floor alpha+1)', flush=True)
