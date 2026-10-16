"""Max number of unmatched demand units n_res = sum_cells (d-c)^+ over class graphs with k E-arcs, subject to
port balance per class and the entropy cut for EVERY union of classes L (expQQ Thm A proof / expQX Step 3).
Matched pairs are neutral for both, so the residual form (min(d,c)=0 per cell) is WLOG; empty classes allowed.

usage: max_nres.py K LAYOUT [MAXM] [ENUM]
  enumerates every placement of K arcs (ordered pairs inside a level, multi-arcs allowed) on LAYOUT,
  maximises n_res; if ENUM=1 also enumerates all residual configs with n_res>=1 and checks residual routing.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, itertools, json
from ortools.sat.python import cp_model
from rg_core import subsets
from route import residual_routable


def model(levels, arcs, maxm, target=None):
    n = len(levels)
    cells = [(a, b) for a in range(n) for b in range(n) if levels[b] == levels[a] + 1]
    md = cp_model.CpModel()
    d = {k: md.NewIntVar(0, maxm, "") for k in cells}
    c = {k: md.NewIntVar(0, maxm, "") for k in cells}
    for k in cells:
        z = md.NewBoolVar("")
        md.Add(d[k] == 0).OnlyEnforceIf(z); md.Add(c[k] == 0).OnlyEnforceIf(z.Not())
    for v in range(n):
        ins = sum(d[k] for k in cells if k[0] == v) + sum(c[k] for k in cells if k[1] == v) + sum(1 for x, y in arcs if y == v)
        outs = sum(d[k] for k in cells if k[1] == v) + sum(c[k] for k in cells if k[0] == v) + sum(1 for x, y in arcs if x == v)
        md.Add(ins == outs)
    for L in subsets(range(n)):
        L = set(L)
        if not L or len(L) == n:
            continue
        lhs = [d[k] for k in cells if k[1] in L and k[0] not in L]
        rhs = [c[k] for k in cells if k[1] in L and k[0] not in L]
        ea = sum(1 for x, y in arcs if y in L and x not in L)
        if lhs:
            md.Add(sum(lhs) <= sum(rhs) + ea)
    obj = sum(d.values())
    if target is None:
        md.Maximize(obj)
    else:
        md.Add(obj >= target)
    return md, d, c, cells


def placements(levels, k):
    n = len(levels)
    pairs = [(x, y) for x in range(n) for y in range(n) if x != y and levels[x] == levels[y]]
    return itertools.combinations_with_replacement(pairs, k)


def main():
    k = int(sys.argv[1]); layout = [int(x) for x in sys.argv[2].split(',')]
    maxm = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    enum = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    levels = []
    for t, m in enumerate(layout):
        levels += [t] * m
    best = (-1, None); hist = {}; nconf = 0; nonroute = []
    for arcs in placements(levels, k):
        md, d, c, cells = model(levels, arcs, maxm)
        s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.parameters.max_time_in_seconds = 60
        r = s.Solve(md)
        assert r == cp_model.OPTIMAL, s.StatusName(r)
        v = int(s.ObjectiveValue())
        hist[v] = hist.get(v, 0) + 1
        if v > best[0]:
            best = (v, arcs, {kk: s.Value(d[kk]) for kk in cells if s.Value(d[kk])}, {kk: s.Value(c[kk]) for kk in cells if s.Value(c[kk])})
        if enum and v >= 1:
            md, d, c, cells = model(levels, arcs, maxm, target=1)
            sols = []

            class CB(cp_model.CpSolverSolutionCallback):
                def on_solution_callback(self):
                    sols.append(({kk: self.Value(d[kk]) for kk in cells if self.Value(d[kk])},
                                 {kk: self.Value(c[kk]) for kk in cells if self.Value(c[kk])}))
            s = cp_model.CpSolver(); s.parameters.enumerate_all_solutions = True; s.parameters.num_workers = 1
            s.Solve(md, CB())
            for dd, cc in sols:
                nconf += 1
                if not residual_routable(levels, arcs, dd, cc):
                    nonroute.append((arcs, dd, cc))
    print(f"k={k} layout={layout} maxm={maxm}: max n_res = {best[0]}  hist={dict(sorted(hist.items()))}  witness={best[1:]}", flush=True)
    if enum:
        print(f"   enumerated {nconf} residual configs with n_res>=1; residual NOT routable: {len(nonroute)}", flush=True)
        for x in nonroute[:5]:
            print("   NONROUTE", x, flush=True)


if __name__ == '__main__':
    main()
