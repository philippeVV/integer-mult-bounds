"""Plane choices r_a, induced orthogonality graph Gamma_E, local search for many edges, CP-SAT alpha."""
import random, sys
from ports import P, dot, rchoices
n = len(P)
NB = [[j for j in range(n) if j != i and dot(P[i], P[j]) == 0] for i in range(n)]
RC = [rchoices(q) for q in P]

def orth(i, ri, j, rj):
    return dot(P[i], P[j]) == 0 and dot(P[i], rj) == 0 and dot(ri, P[j]) == 0 and dot(ri, rj) == 0

def edges(r):
    return [(i, j) for i in range(n) for j in NB[i] if i < j and orth(i, r[i], j, r[j])]

def local_search(seed, iters=40000):
    rng = random.Random(seed)
    r = [rng.choice(RC[i]) for i in range(n)]
    def deg(i, ri): return sum(orth(i, ri, j, r[j]) for j in NB[i])
    cur = sum(deg(i, r[i]) for i in range(n)) // 2
    for it in range(iters):
        i = rng.randrange(n)
        cand = rng.choice(RC[i])
        d = deg(i, cand) - deg(i, r[i])
        if d >= 0 or rng.random() < 0.02:
            r[i] = cand; cur += d
    # final greedy polish
    improved = True
    while improved:
        improved = False
        for i in range(n):
            best = max(RC[i], key=lambda c: deg(i, c))
            if deg(i, best) > deg(i, r[i]):
                r[i] = best; improved = True
    return r

def alpha(r, timeout=60):
    from ortools.sat.python import cp_model
    m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in range(n)]
    for i, j in edges(r): m.AddBoolOr([x[i].Not(), x[j].Not()])
    m.Maximize(sum(x)); s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = timeout
    s.parameters.num_workers = 2
    st = s.Solve(m)
    return int(s.ObjectiveValue()), int(s.BestObjectiveBound()), s.StatusName(st)

if __name__ == '__main__':
    import json
    res = []
    for seed in range(int(sys.argv[1]) if len(sys.argv) > 1 else 4):
        r = local_search(seed)
        e = edges(r)
        a = alpha(r)
        print(seed, 'edges', len(e), 'alpha', a, flush=True)
        res.append(dict(seed=seed, r=r, edges=len(e), alpha=a))
    json.dump(res, open('gamma_runs.json', 'w'))
