"""Rectangle lemma check (two-level class graph).

X(a,b) = #demands a->b  -  #D-arcs a->b  (a in level t, b in level t+1).
Entropy cut with L = (all levels < t) + A' + B' gives, for every A''=A\\A', B':
    X(A'',B') <= cA(A'') + cB(B')
cA(A'') = #E-arcs (inside level t) leaving A'';  cB(B') = #E-arcs (inside t+1) entering B'.
Unmatched (possibly deficient) D-arcs = sum X^- = sum X^+ (totals agree).
Claim tested: max sum X^+ <= kA + kB.  Exhaustive over E-digraphs (multi-arcs allowed), CP-SAT.
"""
import itertools, sys
from ortools.sat.python import cp_model

def cut_out(arcs, S):
    return sum(1 for (x, y) in arcs if x in S and y not in S)

def cut_in(arcs, S):
    return sum(1 for (x, y) in arcs if y in S and x not in S)

def solve(p, q, EA, EB, bound=12):
    md = cp_model.CpModel()
    X = {(a, b): md.NewIntVar(-bound, bound, f"x{a}{b}") for a in range(p) for b in range(q)}
    Xp = {k: md.NewIntVar(0, bound, f"p{k}") for k in X}
    for k in X:
        md.AddMaxEquality(Xp[k], [X[k], 0])
    md.Add(sum(X.values()) == 0)
    for ra in range(1, 1 << p):
        Aset = {a for a in range(p) if ra >> a & 1}
        ca = cut_out(EA, Aset)
        for rb in range(1, 1 << q):
            Bset = {b for b in range(q) if rb >> b & 1}
            md.Add(sum(X[a, b] for a in Aset for b in Bset) <= ca + cut_in(EB, Bset))
    md.Maximize(sum(Xp.values()))
    s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.parameters.max_time_in_seconds = 30
    r = s.Solve(md)
    assert r == cp_model.OPTIMAL, r
    return int(s.ObjectiveValue())

def arcsets(n, k):
    pairs = [(x, y) for x in range(n) for y in range(n) if x != y]
    return itertools.combinations_with_replacement(pairs, k)

worst = 0
for p in range(1, 5):
    for q in range(1, 5):
        if p * q > 12: continue
        for ktot in range(0, 3):   # Lemma R (i)-(ii): f(k) = k proved for k <= 2 (k = 3 reaches 4, see rect_max.py)
            for kA in range(0, ktot + 1):
                kB = ktot - kA
                if (p < 2 and kA) or (q < 2 and kB): continue
                seen = 0
                for EA in arcsets(p, kA):
                    for EB in arcsets(q, kB):
                        v = solve(p, q, EA, EB)
                        seen += 1
                        if v > ktot:
                            print("VIOLATION", p, q, EA, EB, v); sys.exit(1)
                        worst = max(worst, v - ktot)
                print(f"p={p} q={q} kA={kA} kB={kB}: {seen} E-graphs, all max sumX+ <= k", flush=True)
print("max (sumX+ - k) =", worst)
print("PASS_EXPQX_RECT" if worst == 0 else "FAIL_EXPQX_RECT")
