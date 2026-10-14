"""Minimise alpha(Gamma_E) over plane choices (heuristic upper bound on min_r alpha).
Moves: re-plane a vertex of the current MIS to maximise its Gamma_E-degree into the MIS. usage: amin.py seconds seed"""
import sys, json, time, random, numpy as np
from common import *
from ortools.sat.python import cp_model
TL, seed = float(sys.argv[1]), int(sys.argv[2]); rng = random.Random(seed)
r = [rng.choice(RC[a]) for a in range(n)]
def mis(F, hint=None, tl=3.0):
    m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in range(n)]
    ii, jj = np.nonzero(np.triu(F, 1))
    for i, j in zip(ii, jj): m.AddBoolOr([x[i].Not(), x[j].Not()])
    m.Maximize(sum(x)); s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = 1
    s.parameters.random_seed = rng.randrange(10**6); s.Solve(m)
    return [i for i in range(n) if s.Value(x[i])], int(s.BestObjectiveBound())
F = orthmask(r); I, ub = mis(F); best = len(I); t0 = time.time(); it = 0
while time.time() - t0 < TL:
    it += 1
    old = list(r)
    for i in rng.sample(I, min(3, len(I))):
        CM = cand_masks(i, r)
        sc = CM[:, I].sum(1) * 1000 + CM.sum(1) + np.array([rng.random() for _ in range(128)])
        r[i] = RC[i][int(np.argmax(sc))]
    F = orthmask(r); I2, ub2 = mis(F)
    if len(I2) <= len(I) or rng.random() < 0.05: I = I2
    else: r = old; F = orthmask(r)
    if len(I) < best:
        best = len(I); json.dump({'r': [int(x) for x in r], 'alpha_found': best}, open(f'amin_s{seed}.json', 'w'))
        print(f'{time.time()-t0:.0f}s it {it} alpha<= {best} (mis found; edges {F.sum()//2})', flush=True)
print('final', best)
