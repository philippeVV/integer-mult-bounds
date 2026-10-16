"""Part (1)+(2): class graphs for the delta=R hyperbola (k=2 E-arcs, 4 unmatched main climbs).

For each placement of 2 arcs on a layout, find residual configs with n_res = 4 that satisfy port balance and
the RECTANGLE cuts only (expQX Lemma R's hypotheses: L = levels<t  u A'  u B'), then
  (a) list the general node-set entropy cut L that it violates (the obstruction), and
  (b) run exact F2 realisability (one arbitrary-arity gate per class) -- must be UNSAT if (a) is right.
usage: hyper_cg.py LAYOUT [LIMIT]
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, itertools
from ortools.sat.python import cp_model
from rg_core import ClassGraph, realise, subsets
from max_nres import placements


def rect_sets(levels):
    n = len(levels); T = max(levels)
    for t in range(T):
        A = [v for v in range(n) if levels[v] == t]; B = [v for v in range(n) if levels[v] == t + 1]
        low = {v for v in range(n) if levels[v] < t}
        for Ap in subsets(A):
            for Bp in subsets(B):
                yield low | set(Ap) | set(Bp)


def main():
    layout = [int(x) for x in sys.argv[1].split(',')]
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    levels = []
    for t, m in enumerate(layout):
        levels += [t] * m
    n = len(levels)
    cells = [(a, b) for a in range(n) for b in range(n) if levels[b] == levels[a] + 1]
    tot = 0; real = {}
    for arcs in placements(levels, 2):
        md = cp_model.CpModel()
        d = {k: md.NewIntVar(0, 2, "") for k in cells}; c = {k: md.NewIntVar(0, 2, "") for k in cells}
        for k in cells:
            z = md.NewBoolVar(""); md.Add(d[k] == 0).OnlyEnforceIf(z); md.Add(c[k] == 0).OnlyEnforceIf(z.Not())
        for v in range(n):
            ins = sum(d[k] for k in cells if k[0] == v) + sum(c[k] for k in cells if k[1] == v) + sum(1 for x, y in arcs if y == v)
            outs = sum(d[k] for k in cells if k[1] == v) + sum(c[k] for k in cells if k[0] == v) + sum(1 for x, y in arcs if x == v)
            md.Add(ins == outs)
        for L in rect_sets(levels):
            lhs = [d[k] for k in cells if k[1] in L and k[0] not in L]
            if lhs:
                md.Add(sum(lhs) <= sum(c[k] for k in cells if k[1] in L and k[0] not in L)
                       + sum(1 for x, y in arcs if y in L and x not in L))
        md.Add(sum(d.values()) == 4)
        sols = []

        class CB(cp_model.CpSolverSolutionCallback):
            def on_solution_callback(self):
                sols.append(({k: self.Value(d[k]) for k in cells if self.Value(d[k])},
                             {k: self.Value(c[k]) for k in cells if self.Value(c[k])}))
                if len(sols) >= limit:
                    self.StopSearch()
        s = cp_model.CpSolver(); s.parameters.enumerate_all_solutions = True; s.parameters.num_workers = 1
        s.Solve(md, CB())
        for dd, cc in sols:
            tot += 1
            G = ClassGraph(levels, arcs, dd, cc)
            viol = None
            for L in subsets(range(n)):
                L = set(L)
                lhs = sum(m for (a, b), m in dd.items() if b in L and a not in L)
                rhs = sum(m for (a, b), m in cc.items() if b in L and a not in L) + sum(1 for x, y in arcs if y in L and x not in L)
                if lhs > rhs:
                    viol = (sorted(L), lhs, rhs); break
            r, _ = realise(G, 2, time_limit=60)
            real[r] = real.get(r, 0) + 1
            if tot <= 6:
                print(f"arcs={arcs} demands={dd} climbs={cc}\n   violated cut L={viol}  F2: {r}", flush=True)
    print(f"layout={layout}: {tot} rectangle-feasible n_res=4 class graphs; F2 realisability {real}", flush=True)
    print("PASS_EXPRG_HYPER" if set(real) <= {"UNSAT"} else "FAIL_EXPRG_HYPER")   # no realisable n_res = 4 class graph with 2 arcs


if __name__ == '__main__':
    main()
