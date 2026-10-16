"""expRG core: class graphs with levels, entropy cuts, and exact F2/Fp realisability by CP-SAT.

A class graph has classes (cid -> level), E-arcs (tail, head) inside a level, and per adjacent level pair
multiplicities d[a,b] (demands: a source in class a whose sink is in class b) and c[a,b] (main climbs a->b).
Each class is modelled as ONE gate of arbitrary arity (most permissive F-realisation of a convex class).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import itertools
from ortools.sat.python import cp_model


def subsets(xs):
    xs = list(xs)
    for r in range(len(xs) + 1):
        yield from itertools.combinations(xs, r)


class ClassGraph:
    def __init__(self, levels, arcs, d, c):
        self.levels = levels            # list: class -> level
        self.arcs = list(arcs)          # (tail, head) same level
        self.d = dict(d)                # (a,b) -> count, level(b)=level(a)+1
        self.c = dict(c)

    def ports(self, v):
        src = sum(n for (a, b), n in self.d.items() if a == v)
        snk = sum(n for (a, b), n in self.d.items() if b == v)
        cin = sum(n for (a, b), n in self.c.items() if b == v)
        cout = sum(n for (a, b), n in self.c.items() if a == v)
        ein = sum(1 for (x, y) in self.arcs if y == v)
        eout = sum(1 for (x, y) in self.arcs if x == v)
        return src + cin + ein, snk + cout + eout

    def balanced(self):
        return all(i == o and i > 0 for i, o in (self.ports(v) for v in range(len(self.levels))))

    def cut_ok(self):
        n = len(self.levels)
        for L in subsets(range(n)):
            L = set(L)
            lhs = sum(m for (a, b), m in self.d.items() if b in L and a not in L)
            rhs = sum(m for (a, b), m in self.c.items() if b in L and a not in L)
            rhs += sum(1 for (x, y) in self.arcs if y in L and x not in L)
            if lhs > rhs:
                return False
        return True

    def nres(self):
        keys = set(self.d) | set(self.c)
        return sum(max(self.d.get(k, 0) - self.c.get(k, 0), 0) for k in keys)


def realise(G, p=2, time_limit=60.0, workers=1):
    """Exact existence over F_p (p=2 via XOR, p prime small via modular ints) of gate maps such that the
    class-gate circuit computes a permutation consistent with the class-level demands.
    Wire values are linear forms in the W source variables. Gate outputs = M * gate inputs (any M; invertibility
    is forced because the sinks must be W distinct unit vectors). Returns (status, solver)."""
    n = len(G.levels)
    # sources
    srcs = []  # (class, sink-class)
    for (a, b), m in sorted(G.d.items()):
        srcs += [(a, b)] * m
    W = len(srcs)
    md = cp_model.CpModel()
    # wires: inputs to each class, outputs of each class
    inw = {v: [] for v in range(n)}
    outw = {v: [] for v in range(n)}
    # source wires: fixed unit forms
    for w, (a, b) in enumerate(srcs):
        inw[a].append(('const', [1 if j == w else 0 for j in range(W)]))
    # topological order: by level, then arcs DAG inside level
    order = topo(G)
    splits = []
    if order is None:
        # antiparallel pairs x->y, y->x inside a level: split x into x1 (own inputs + other incoming arcs ->
        # arc to y + j internal wires) and x2 (internal wires + arc back from y -> all exits + other outgoing
        # arcs), j = #inputs(x1) - 1. Canonical (most-information) split; the result must be acyclic.
        seen = set()
        for (x, y) in G.arcs:
            if (y, x) in G.arcs and x not in seen and y not in seen:
                splits.append((x, y)); seen |= {x, y}
        if not splits:
            return 'CYCLIC', None
    split_of = {x: y for x, y in splits}
    # create edge variables for climbs and arcs (form vectors)
    def newform(tag):
        if p == 2:
            return [md.NewBoolVar(f"{tag}_{j}") for j in range(W)]
        return [md.NewIntVar(0, p - 1, f"{tag}_{j}") for j in range(W)]
    k = 0
    for (a, b), m in sorted(G.c.items()):
        for _ in range(m):
            f = newform(f"c{k}"); k += 1
            outw[a].append(('var', f)); inw[b].append(('var', f))
    arcf = {}; back = {}
    gedges = []   # gate-level DAG edges for the acyclicity check
    for (x, y) in G.arcs:
        f = newform(f"e{k}"); k += 1
        if y in split_of and split_of[y] == x:            # back arc into split class y -> y2
            outw[x].append(('var', f)); back.setdefault(y, []).append(('var', f)); continue
        if x in split_of and split_of[x] == y and x not in arcf:  # forward arc from split class -> from x1
            arcf[x] = ('var', f); inw[y].append(arcf[x]); continue
        outw[x].append(('var', f)); inw[y].append(('var', f))
    for w, (a, b) in enumerate(srcs):
        outw[b].append(('const', [1 if j == w else 0 for j in range(W)]))
    gates = list(range(n))
    for x, y in splits:
        j = len(inw[x]) - 1
        mids = [('var', newform(f"m{x}_{i}")) for i in range(j)]
        x2 = ('x2', x)
        outw[x2] = outw[x]
        outw[x] = [arcf[x]] + mids
        inw[x2] = mids + back.get(x, [])
        gates.append(x2)
    # acyclicity at gate level (wires shared by identity of form lists)
    prod = {}
    for v in gates:
        for (t, f) in outw[v]:
            if t == 'var':
                prod[id(f)] = v
    succ = {v: set() for v in gates}
    for v in gates:
        for (t, f) in inw[v]:
            if t == 'var':
                succ[prod[id(f)]].add(v)
    indeg = {v: 0 for v in gates}
    for v in gates:
        for u in succ[v]:
            indeg[u] += 1
    ready = [v for v in gates if indeg[v] == 0]; cnt = 0
    while ready:
        v = ready.pop(); cnt += 1
        for u in succ[v]:
            indeg[u] -= 1
            if indeg[u] == 0:
                ready.append(u)
    if cnt != len(gates):
        return 'CYCLIC', None
    # gate constraints: each output = sum_i coef_i * input_i
    for v in gates:
        ins, outs = inw[v], outw[v]
        assert len(ins) == len(outs), (v, len(ins), len(outs))
        for oi, (ot, of) in enumerate(outs):
            if p == 2:
                coef = [md.NewBoolVar(f"g{v}_{oi}_{ii}") for ii in range(len(ins))]
            else:
                coef = [md.NewIntVar(0, p - 1, f"g{v}_{oi}_{ii}") for ii in range(len(ins))]
            for j in range(W):
                terms = []
                for ii, (it, iform) in enumerate(ins):
                    x = iform[j]
                    if it == 'const':
                        if x:
                            terms.append(coef[ii])
                    else:
                        if p == 2:
                            t = md.NewBoolVar("")
                            md.AddBoolAnd([coef[ii], x]).OnlyEnforceIf(t)
                            md.AddBoolOr([coef[ii].Not(), x.Not()]).OnlyEnforceIf(t.Not())
                            terms.append(t)
                        else:
                            t = md.NewIntVar(0, (p - 1) ** 2, "")
                            md.AddMultiplicationEquality(t, [coef[ii], x])
                            terms.append(t)
                tgt = of[j]
                if p == 2:
                    if ot == 'const':
                        if terms:
                            md.AddBoolXOr(terms + ([md.NewConstant(1)] if tgt == 0 else []))
                        elif tgt:
                            return 'UNSAT', None
                    else:
                        md.AddBoolXOr(terms + [tgt.Not()] if terms else [tgt.Not()])
                else:
                    q = md.NewIntVar(0, 10 * len(ins) * p * p, "")
                    rhs = tgt if ot == 'var' else int(tgt)
                    md.Add(sum(terms) == q * p + rhs)
    s = cp_model.CpSolver()
    s.parameters.num_workers = workers
    s.parameters.max_time_in_seconds = time_limit
    r = s.Solve(md)
    return {cp_model.OPTIMAL: 'SAT', cp_model.FEASIBLE: 'SAT', cp_model.INFEASIBLE: 'UNSAT'}.get(r, 'UNKNOWN'), s


def topo(G):
    n = len(G.levels)
    succ = {v: set() for v in range(n)}
    for (x, y) in G.arcs:
        succ[x].add(y)
    indeg = {v: 0 for v in range(n)}
    for x in succ:
        for y in succ[x]:
            indeg[y] += 1
    out = []
    ready = sorted([v for v in range(n) if indeg[v] == 0], key=lambda v: G.levels[v])
    while ready:
        v = ready.pop(0); out.append(v)
        for y in succ[v]:
            indeg[y] -= 1
            if indeg[y] == 0:
                ready.append(y)
        ready.sort(key=lambda v: G.levels[v])
    return out if len(out) == n else None
