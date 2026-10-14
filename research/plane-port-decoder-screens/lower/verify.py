"""Exact checks for the plane-port decoder floor notes (MODEL.md). Standard library + optional ortools
(only to *find* independent sets; every floor uses an exactly re-checked witness set)."""
import json, itertools, sys
from ports import P, dot, rchoices
from gamma import edges, NB, n, orth
ok = lambda c, m: print(('PASS ' if c else 'FAIL ') + m) or (c or sys.exit(1))
# 1. ports
ok(len(P) == 112 and all(dot(q, q) for q in P) and {bin(q).count('1') for q in P} == {3, 7}, 'ports: 112 odd vectors of weight 3/7')
ok(all(len(NB[i]) == 52 for i in range(n)) and sum(map(len, NB)) == 5824, 'P112 graph 52-regular, 2912 edges')
miss = sorted(x for x in range(512) if bin(x).count('1') in (3, 7) and x not in P)
ok(len(miss) == 8 and all(dot(a, b) == 1 for a, b in itertools.combinations(miss, 2)), 'missing 8 = frame (pairwise dot 1)')
# 2. plane normal form: nondegenerate planes through q <-> r in q^perp with r.r = 1, unique such r per plane
for q in P:
    planes = {}
    for r in range(512):
        if r in (0, q): continue
        pl = frozenset((0, q, r, q ^ r))
        g = [[dot(a, b) for b in (q, r)] for a in (q, r)]
        if (g[0][0] * g[1][1] + g[0][1] * g[1][0]) % 2: planes.setdefault(pl, []).append(r)
    reps = [[r for r in rs if dot(r, q) == 0 and dot(r, r) == 1] for rs in planes.values()]
    assert len(planes) == 128 and all(len(x) == 1 for x in reps) and sorted(x[0] for x in reps) == rchoices(q)
ok(True, 'each port: exactly 128 nondegenerate planes, each with unique orthonormal basis {q, r}')
# 3. Gamma_E subset of P112 and degree bound 28
mx = [max(sum(1 for j in NB[i] if dot(P[j], r) == 0) for r in rchoices(P[i])) for i in range(n)]
ok(set(mx) == {28}, 'every port: at most 28 of its 52 P112 neighbours can stay orthogonal (any planes); |Gamma_E| <= 1568')
# 4. instances: witness independent sets give the rigorous floor C >= |S| + 1
def witness(r):
    try:
        from ortools.sat.python import cp_model
    except ImportError:
        return None
    m = cp_model.CpModel(); x = [m.NewBoolVar('') for _ in range(n)]
    for i, j in edges(r): m.AddBoolOr([x[i].Not(), x[j].Not()])
    m.Maximize(sum(x)); s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = 60; s.parameters.num_workers = 2
    st = s.Solve(m)
    return [i for i in range(n) if s.Value(x[i])], s.StatusName(st), int(s.BestObjectiveBound())
for f in sys.argv[1:]:
    d = json.load(open(f)); runs = d if isinstance(d, list) else [d]
    for run in runs[:1]:
        r = run['r']; E = set(edges(r))
        assert all(dot(r[i], P[i]) == 0 and dot(r[i], r[i]) == 1 for i in range(n))
        w = witness(r)
        if w is None: print('SKIP witness (no ortools)'); continue
        S, st, ub = w
        ok(all((min(a, b), max(a, b)) not in E for a, b in itertools.combinations(S, 2)),
           f'{f}: |Gamma_E|={len(E)}; witness set of {len(S)} pairwise non-orthogonal planes (solver {st}, alpha<= {ub}) => floor C >= {len(S)+1}')
print('DONE')
