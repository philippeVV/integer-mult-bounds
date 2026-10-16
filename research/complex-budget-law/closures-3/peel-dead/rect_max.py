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

worst = 0; top = {}
for p in range(1, 5):
    for q in range(1, 5):
        for ktot in range(0, 4):
            if p * q > 16 or (p * q > 9 and ktot > 2): continue
            for kA in range(0, ktot + 1):
                kB = ktot - kA
                if (p < 2 and kA) or (q < 2 and kB): continue
                seen = 0; best=[-1,None,None]
                for EA in arcsets(p, kA):
                    for EB in arcsets(q, kB):
                        v = solve(p, q, EA, EB)
                        seen += 1
                        if v > best[0]: best[:] = [v, EA, EB]
                worst = max(worst, best[0] - ktot); top[ktot] = max(top.get(ktot, -1), best[0])
                print(f"p={p} q={q} kA={kA} kB={kB}: {seen} graphs max sumX+ = {best[0]} (k={ktot}) at {best[1]} {best[2]}", flush=True)
print("max (sumX+ - k) =", worst, "; max sumX+ by k:", top)
# Lemma R: f(k) = k for k <= 2; general bound kA + kB + 2 kA kB.  Expected maxima (expQX): k = 0, 1, 2 -> k; k = 3 -> 4
print("PASS_EXPQX_RECTMAX" if top == {0: 0, 1: 1, 2: 2, 3: 4} else "FAIL_EXPQX_RECTMAX")
