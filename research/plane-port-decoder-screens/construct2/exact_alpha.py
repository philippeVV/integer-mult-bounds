import sys, json, numpy as np
from common import *
from ortools.sat.python import cp_model
for f in sys.argv[1:]:
    r = json.load(open(f))['r']; F = orthmask(r)
    m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in range(n)]
    for i, j in zip(*np.nonzero(np.triu(F, 1))): m.AddBoolOr([x[i].Not(), x[j].Not()])
    m.Maximize(sum(x)); s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = 240; s.parameters.num_workers = 2
    st = s.Solve(m); print(f, 'alpha', int(s.ObjectiveValue()), 'bound', int(s.BestObjectiveBound()), s.StatusName(st), 'edges', F.sum() // 2, flush=True)
