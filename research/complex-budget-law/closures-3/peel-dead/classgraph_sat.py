"""CP-SAT violation search over graded class graphs (expQX §1 model).

Classes (t,i): level t in 0..T-1, i in 0..p-1.  Arcs:
  D-arcs   (main climbs)  a(t) -> b(t+1), multiplicity d[a,b] in 0..MU
  demands                 a(t) -> b(t+1), multiplicity q[a,b] in 0..MU
  E_j  (j<k, each placed exactly once) u(t) -> v(t), u != v
Potential c[v] in Z^k (E-coordinates of the frame, M_v = M_0 + l(v) I + c(v).E):
  E_j arc u->v  =>  c[v]-c[u] = e_j ;  demand a->b  =>  c[b] = c[a] ;  D a->b has P_D = c[b]-c[a].
(I, E_1..E_k Q-independent; the potential is exactly the frame data, so this is no restriction.)

Hypotheses (flag --cuts):
  rect : level equalities  sum dem(t->t+1) = sum D(t->t+1)  + all rectangle cuts
         X(A'',B') <= cE_out_t(A'') + cE_in_{t+1}(B')   for all A'' in level t, B' in level t+1
         (exactly the hypotheses used in the proof; weaker than the full entropy cut)
  all  : the entropy cut for EVERY node set L:  #dem into L from outside <= #climbs entering L
  none : only level equalities (negative control)
Optional --circ: route conservation  (climb in - climb out = dem in - dem out at every class).

Claims (--claim), each searched as a violation (INFEASIBLE = claim holds in the box):
  far      : some D a->b at pair t has P_D[j] != 0 with level(E_j) outside {t-1..t+2}  (the stop test)
  far0     : ... with level(E_j) outside {t, t+1}
  comp     : ... with E_j outside comp(a) u comp(b) (E-components), certified by E-closed sets
  balance  : X(K,K') != 0 for some union of E-components K (level t), K' (level t+1)
  use2     : some E_j has P_D[j] != 0 for >= 2 D-arcs (with multiplicity) of one adjacent pair
  use3     : some E_j has P_D[j] != 0 for >= 3 D-arcs in total
  coef     : some D has |P_D[j]| >= 2
"""
import argparse, itertools, sys, time
from ortools.sat.python import cp_model

ap = argparse.ArgumentParser()
ap.add_argument('--T', type=int, default=5)
ap.add_argument('--p', type=int, default=3)
ap.add_argument('--k', type=int, default=4)
ap.add_argument('--mu', type=int, default=3)
ap.add_argument('--cb', type=int, default=6)
ap.add_argument('--cuts', default='rect')
ap.add_argument('--circ', action='store_true')
ap.add_argument('--claim', default='far')
ap.add_argument('--time', type=float, default=600)
ap.add_argument('--workers', type=int, default=1)
ap.add_argument("--thr", type=int, default=3)
ap.add_argument("--expect", default="INFEASIBLE")   # closures-3: PASS marker when the solver status matches
a = ap.parse_args()
T, p, k, MU, CB = a.T, a.p, a.k, a.mu, a.cb

md = cp_model.CpModel()
V = [(t, i) for t in range(T) for i in range(p)]
P2 = [(t, a_, b_) for t in range(T - 1) for a_ in range(p) for b_ in range(p)]  # (t, a in t, b in t+1)
d = {x: md.NewIntVar(0, MU, f"d{x}") for x in P2}
q = {x: md.NewIntVar(0, MU, f"q{x}") for x in P2}
hd = {x: md.NewBoolVar(f"hd{x}") for x in P2}
hq = {x: md.NewBoolVar(f"hq{x}") for x in P2}
for x in P2:
    md.Add(d[x] >= 1).OnlyEnforceIf(hd[x]); md.Add(d[x] == 0).OnlyEnforceIf(hd[x].Not())
    md.Add(q[x] >= 1).OnlyEnforceIf(hq[x]); md.Add(q[x] == 0).OnlyEnforceIf(hq[x].Not())
ES = [(t, u, v) for t in range(T) for u in range(p) for v in range(p) if u != v]
y = {(j,) + s: md.NewBoolVar(f"y{j}{s}") for j in range(k) for s in ES}
lev = {}
for j in range(k):
    md.AddExactlyOne(y[(j,) + s] for s in ES)
    lev[j] = md.NewIntVar(0, T - 1, f"lev{j}")
    md.Add(lev[j] == sum(s[0] * y[(j,) + s] for s in ES))
c = {(v, j): md.NewIntVar(-CB, CB, f"c{v}{j}") for v in V for j in range(k)}
for j in range(k):
    for (t, u, v) in ES:
        for jj in range(k):
            md.Add(c[(t, v), jj] - c[(t, u), jj] == (1 if jj == j else 0)).OnlyEnforceIf(y[(j, t, u, v)])
for (t, a_, b_) in P2:
    for jj in range(k):
        md.Add(c[(t + 1, b_), jj] == c[(t, a_), jj]).OnlyEnforceIf(hq[(t, a_, b_)])
md.Add(sum(d.values()) >= 1)

def Eout(t, S):  # E-arcs at level t leaving S (S subset of level-t indices)
    return sum(y[(j, t, u, v)] for j in range(k) for u in S for v in range(p) if v not in S and u != v)
def Ein(t, S):
    return sum(y[(j, t, u, v)] for j in range(k) for v in S for u in range(p) if u not in S and u != v)

for t in range(T - 1):
    md.Add(sum(q[(t, a_, b_)] for a_ in range(p) for b_ in range(p)) ==
           sum(d[(t, a_, b_)] for a_ in range(p) for b_ in range(p)))
subsets = [set(s) for r in range(1, p + 1) for s in itertools.combinations(range(p), r)]
def cut_L(L):  # entropy cut for node set L (only arcs touching L matter)
    lhs = sum(q[(t, a_, b_)] for (t, a_, b_) in P2 if (t + 1, b_) in L and (t, a_) not in L)
    rhs = sum(d[(t, a_, b_)] for (t, a_, b_) in P2 if (t + 1, b_) in L and (t, a_) not in L)
    rhs += sum(y[(j, t, u, v)] for j in range(k) for (t, u, v) in ES if (t, v) in L and (t, u) not in L)
    md.Add(lhs <= rhs)
if a.cuts == 'stair3':
    # three-level staircase cuts L = (levels < t-1) u S_{t-1} u S_t u S_{t+1}  (levels > t+1 empty);
    # these include the rectangle cuts (S_{t-1} = whole level) and the level cuts
    allsub = [set(s) for r in range(0, p + 1) for s in itertools.combinations(range(p), r)]
    for t in range(1, T - 1):
        for S0 in allsub:
            for S1 in allsub:
                for S2 in allsub:
                    L = {(tt, i) for tt in range(t - 1) for i in range(p)}
                    L |= {(t - 1, i) for i in S0} | {(t, i) for i in S1} | {(t + 1, i) for i in S2}
                    cut_L(L)
    for S1 in allsub:      # bottom pair (t-1 = -1)
        for S2 in allsub:
            cut_L({(0, i) for i in S1} | {(1, i) for i in S2})
    if a.cuts == 'stair3': pass
if a.cuts == 'rect':
    for t in range(T - 1):
        for A in subsets:
            for B in subsets:
                md.Add(sum(q[(t, x, z)] - d[(t, x, z)] for x in A for z in B) <= Eout(t, A) + Ein(t + 1, B))
elif a.cuts == 'all':
    n = len(V)
    for mask in range(1, 1 << n):
        L = {V[i] for i in range(n) if mask >> i & 1}
        lhs = sum(q[(t, a_, b_)] for (t, a_, b_) in P2 if (t + 1, b_) in L and (t, a_) not in L)
        rhs = sum(d[(t, a_, b_)] for (t, a_, b_) in P2 if (t + 1, b_) in L and (t, a_) not in L)
        rhs += sum(y[(j, t, u, v)] for j in range(k) for (t, u, v) in ES if (t, v) in L and (t, u) not in L)
        md.Add(lhs <= rhs)
if a.circ:
    for (t, i) in V:
        din = sum(d[(t - 1, x, i)] for x in range(p)) if t > 0 else 0
        dout = sum(d[(t, i, x)] for x in range(p)) if t < T - 1 else 0
        qin = sum(q[(t - 1, x, i)] for x in range(p)) if t > 0 else 0
        qout = sum(q[(t, i, x)] for x in range(p)) if t < T - 1 else 0
        ein = sum(y[(j, t, u, i)] for j in range(k) for u in range(p) if u != i)
        eout = sum(y[(j, t, i, v)] for j in range(k) for v in range(p) if v != i)
        md.Add(din + ein - dout - eout == qin - qout)

def ne_bool(x, j):  # bool: P_D[j] != 0 for D-slot x
    (t, a_, b_) = x
    z = md.NewBoolVar("")
    diff = c[(t + 1, b_), j] - c[(t, a_), j]
    md.Add(diff != 0).OnlyEnforceIf(z); md.Add(diff == 0).OnlyEnforceIf(z.Not())
    return z

def closed_set(t, name):  # bools on level t, closed under E-arcs (union of E-components)
    s = [md.NewBoolVar(f"{name}{i}") for i in range(p)]
    for (tt, u, v) in ES:
        if tt != t: continue
        for j in range(k):
            md.AddImplication(y[(j, t, u, v)], s[u]).OnlyEnforceIf(s[v])
            md.AddImplication(y[(j, t, u, v)], s[v]).OnlyEnforceIf(s[u])
    return s

cl = a.claim
if cl in ('far', 'far0', 'comp', 'coef'):
    t = md.NewIntVar(0, T - 2, "t"); j_ = md.NewIntVar(0, k - 1, "j")
    pick = {x: md.NewBoolVar(f"pk{x}") for x in P2}
    md.AddExactlyOne(pick.values())
    pj = [md.NewBoolVar(f"pj{j}") for j in range(k)]
    md.AddExactlyOne(pj)
    for x in P2:
        md.AddImplication(pick[x], hd[x])
        md.Add(t == x[0]).OnlyEnforceIf(pick[x])
    for x in P2:
        for j in range(k):
            both = md.NewBoolVar("")
            md.AddBoolAnd([pick[x], pj[j]]).OnlyEnforceIf(both)
            md.AddBoolOr([pick[x].Not(), pj[j].Not()]).OnlyEnforceIf(both.Not())
            diff = c[(x[0] + 1, x[2]), j] - c[(x[0], x[1]), j]
            if cl == 'coef':
                md.AddAbsEquality(md.NewIntVar(2, 2 * CB, ""), diff).OnlyEnforceIf(both)
            else:
                md.Add(diff != 0).OnlyEnforceIf(both)
            if cl == 'far':
                fl = md.NewBoolVar("")
                md.Add(lev[j] <= x[0] - 2).OnlyEnforceIf([both, fl])
                md.Add(lev[j] >= x[0] + 3).OnlyEnforceIf([both, fl.Not()])
            elif cl == 'far0':
                fl = md.NewBoolVar("")
                md.Add(lev[j] <= x[0] - 1).OnlyEnforceIf([both, fl])
                md.Add(lev[j] >= x[0] + 2).OnlyEnforceIf([both, fl.Not()])
            elif cl == 'comp':
                # certify E_j not in comp(a) and not in comp(b): E-closed sets Sa (level t) containing a,
                # Sb (level t+1) containing b, with E_j's arc outside both
                pass
    if cl == 'comp':
        Sa = {tt: closed_set(tt, f"Sa{tt}_") for tt in range(T - 1)}
        Sb = {tt: closed_set(tt + 1, f"Sb{tt}_") for tt in range(T - 1)}
        for x in P2:
            (tt, a_, b_) = x
            md.AddImplication(pick[x], Sa[tt][a_]); md.AddImplication(pick[x], Sb[tt][b_])
            for j in range(k):
                for (lt, u, v) in ES:
                    # if E_j sits at level tt inside Sa, or level tt+1 inside Sb -> forbidden
                    if lt == tt:
                        md.AddBoolOr([pick[x].Not(), pj[j].Not(), y[(j, lt, u, v)].Not(), Sa[tt][u].Not()])
                    if lt == tt + 1:
                        md.AddBoolOr([pick[x].Not(), pj[j].Not(), y[(j, lt, u, v)].Not(), Sb[tt][u].Not()])
elif cl == 'balance':
    tsel = [md.NewBoolVar(f"ts{t}") for t in range(T - 1)]
    md.AddExactlyOne(tsel)
    S = {t: (closed_set(t, f"A{t}_"), closed_set(t + 1, f"B{t}_")) for t in range(T - 1)}
    for t in range(T - 1):
        A, B = S[t]
        terms = []
        for a_ in range(p):
            for b_ in range(p):
                w = md.NewBoolVar("")
                md.AddBoolAnd([A[a_], B[b_]]).OnlyEnforceIf(w)
                md.AddBoolOr([A[a_].Not(), B[b_].Not()]).OnlyEnforceIf(w.Not())
                xv = md.NewIntVar(-MU, MU, "")
                md.Add(xv == q[(t, a_, b_)] - d[(t, a_, b_)]).OnlyEnforceIf(w)
                md.Add(xv == 0).OnlyEnforceIf(w.Not())
                terms.append(xv)
        neq = md.NewBoolVar("")
        md.Add(sum(terms) != 0).OnlyEnforceIf(tsel[t])
elif cl in ('use2', 'use3'):
    jsel = [md.NewBoolVar(f"js{j}") for j in range(k)]
    md.AddExactlyOne(jsel)
    tot = []
    per_t = {}
    for x in P2:
        for j in range(k):
            z = ne_bool(x, j)
            w = md.NewIntVar(0, MU, "")
            md.Add(w == d[x]).OnlyEnforceIf([z, jsel[j]])
            md.Add(w == 0).OnlyEnforceIf(z.Not()); md.Add(w == 0).OnlyEnforceIf(jsel[j].Not())
            per_t.setdefault(x[0], []).append(w); tot.append(w)
    if cl == 'use2':
        tsel = [md.NewBoolVar(f"ts{t}") for t in range(T - 1)]
        md.AddExactlyOne(tsel)
        for t in range(T - 1):
            md.Add(sum(per_t[t]) >= 2).OnlyEnforceIf(tsel[t])
    else:
        md.Add(sum(tot) >= a.thr)

s = cp_model.CpSolver(); s.parameters.num_workers = a.workers; s.parameters.max_time_in_seconds = a.time
t0 = time.time(); r = s.Solve(md)
print(f"T={T} p={p} k={k} mu={MU} cb={CB} cuts={a.cuts} circ={a.circ} claim={cl}: {s.StatusName(r)} ({time.time()-t0:.1f}s)", flush=True)
if r in (cp_model.FEASIBLE, cp_model.OPTIMAL):
    for j in range(k):
        for sv in ES:
            if s.Value(y[(j,) + sv]): print(f"  E{j}: level {sv[0]}  {sv[1]}->{sv[2]}")
    for x in P2:
        if s.Value(d[x]) or s.Value(q[x]):
            P = [s.Value(c[(x[0] + 1, x[2]), j]) - s.Value(c[(x[0], x[1]), j]) for j in range(k)]
            print(f"  pair {x[0]}: {x[1]}->{x[2]}  D x{s.Value(d[x])} dem x{s.Value(q[x])}  P_D={P}")
st_ = 'FEASIBLE' if s.StatusName(r) in ('FEASIBLE', 'OPTIMAL') else s.StatusName(r)   # OPTIMAL also means a violation was found
print(('PASS_EXPRH_' if st_ == a.expect else 'FAIL_EXPRH_') + cl.upper() + '_' + a.expect)
