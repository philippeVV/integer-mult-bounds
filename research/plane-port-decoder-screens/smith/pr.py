from common import *
from ortools.sat.python import cp_model
def prk(Ea, Eb): return rank2([sum(dot(x, y) << j for j, y in enumerate(Eb)) for x in Ea])
def PRmat(r):
    E = planes(r); return [[prk(E[a], E[b]) for b in range(n)] for a in range(n)]
def mis(adj_ok, t=60, workers=2):
    """max set with all pairs satisfying adj_ok(a,b) (i.e. independent in complement)."""
    m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            if not adj_ok(a, b): m.AddBoolOr([x[a].Not(), x[b].Not()])
    m.Maximize(sum(x)); s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = t; s.parameters.num_workers = workers
    st = s.Solve(m)
    return [a for a in range(n) if s.Value(x[a])], int(s.BestObjectiveBound()), s.StatusName(st)
