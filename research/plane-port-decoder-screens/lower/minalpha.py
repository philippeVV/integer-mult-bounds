"""Local search minimising alpha(Gamma_E) over plane choices (heuristic; upper bound on min alpha)."""
import random, sys, json
from gamma import *
from ortools.sat.python import cp_model

def mis(r, timeout=10):
    m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in range(n)]
    for i, j in edges(r): m.AddBoolOr([x[i].Not(), x[j].Not()])
    m.Maximize(sum(x)); s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = timeout; s.parameters.num_workers = 2
    s.Solve(m)
    return [i for i in range(n) if s.Value(x[i])], int(s.BestObjectiveBound())

seed = int(sys.argv[1]); rng = random.Random(seed)
r = [rng.choice(RC[i]) for i in range(n)]
S, ub = mis(r); best = len(S); bestr = list(r)
for it in range(int(sys.argv[2])):
    # pick a vertex in the MIS and re-choose its plane to maximise orth neighbours inside S
    i = rng.choice(S)
    Sset = set(S)
    scores = [(sum(orth(i, c, j, r[j]) for j in NB[i] if j in Sset), sum(orth(i, c, j, r[j]) for j in NB[i]), c) for c in RC[i]]
    mx = max(s[:2] for s in scores)
    cands = [c for a, b, c in scores if (a, b) == mx]
    old = r[i]; r[i] = rng.choice(cands)
    S2, ub2 = mis(r)
    if len(S2) <= len(S) or rng.random() < 0.1:
        S = S2
    else:
        r[i] = old
    if len(S) < best:
        best = len(S); bestr = list(r); json.dump(dict(seed=seed, alpha=best, r=bestr), open(f"minalpha_{seed}.json", "w")); print(it, 'alpha', best, 'edges', len(edges(r)), flush=True)
json.dump(dict(seed=seed, alpha=best, r=bestr), open(f'minalpha_{seed}.json', 'w'))
print('final', best)
