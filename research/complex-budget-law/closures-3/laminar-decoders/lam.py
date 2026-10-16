"""expQS lam.py: laminar (helper-free) order of direct corrections, OR-reach frame semantics.

Registers: source x_b (starts at frame <p_b> plus optional extra labels), target t_a (starts at optional labels).
A correction edge e = (b, a) is the direct gate t_a += c * x_b; it needs (a, b) free (B[a][b] = 0).
A gate joins frames: after it both wires carry the union of the label sets they carried before (OR-reach).
Target t_a's frame must stay inside p_a^perp at all times: every label c that ever reaches t_a, or that x_b
carries when it is read into t_a, must satisfy B[c][a] = 0 (c = a is always forbidden: B[a][a] = 1).
Laminar order = a total order of the edges for which no forbidden label arrives (the reader-chain lemma of
expNT generalised to arbitrary correction graphs: T_{a_i} perp U_{a_k} for earlier readers a_i of a source).

check(B, edges, order)   exact replay of an order (returns first violation or None)
sat(B, edges, ...)       CP-SAT: exists a laminar order?  (integer times + monotone reach booleans)
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import itertools


def replay(B, edges, order, x0=None, t0=None):
    xs = {}; ts = {}
    for b, a in edges:
        xs.setdefault(b, {b} | set((x0 or {}).get(b, ())))
        ts.setdefault(a, set((t0 or {}).get(a, ())))
    for a, S in ts.items():
        for c in S:
            if B[c][a]:
                return ('t0', a, c)
    for i in order:
        b, a = edges[i]
        U = xs[b] | ts[a]
        bad = [c for c in U if B[c][a]]
        if bad:
            return (i, edges[i], bad[0])
        xs[b] = U; ts[a] = set(U)
    return None


def brute(B, edges, x0=None, t0=None, limit=None):
    for n, order in enumerate(itertools.permutations(range(len(edges)))):
        if limit and n > limit:
            return 'unknown'
        if replay(B, edges, order, x0, t0) is None:
            return order
    return None


def sat(B, edges, x0=None, t0=None, time=60, workers=2, labels=None, hint=None):
    from ortools.sat.python import cp_model
    E = len(edges)
    m = cp_model.CpModel()
    T = [m.NewIntVar(0, E - 1, f't{i}') for i in range(E)]
    m.AddAllDifferent(T)
    # adjacency: edges sharing source wire, edges sharing target wire
    bysrc = {}; bytgt = {}
    for i, (b, a) in enumerate(edges):
        bysrc.setdefault(b, []).append(i); bytgt.setdefault(a, []).append(i)
    before = {}
    def bef(f, e):
        if (f, e) not in before:
            v = m.NewBoolVar(f'o{f}_{e}')
            m.Add(T[f] < T[e]).OnlyEnforceIf(v); m.Add(T[f] > T[e]).OnlyEnforceIf(v.Not())
            before[(f, e)] = v; before[(e, f)] = v.Not()
        return before[(f, e)]
    ports = sorted({b for b, a in edges} | {a for b, a in edges})
    labs = labels if labels is not None else ports
    # carried label sets just before edge e on its source wire (Rx) and target wire (Rt), at-least semantics
    Rx = {}; Rt = {}
    for i, (b, a) in enumerate(edges):
        for c in labs:
            Rx[c, i] = m.NewBoolVar(''); Rt[c, i] = m.NewBoolVar('')
            init_x = c == b or c in (x0 or {}).get(b, ())
            init_t = c in (t0 or {}).get(a, ())
            if init_x: m.Add(Rx[c, i] == 1)
            if init_t: m.Add(Rt[c, i] == 1)
            if B[c][a]:
                m.Add(Rx[c, i] == 0); m.Add(Rt[c, i] == 0)
    for i, (b, a) in enumerate(edges):
        for f in bysrc[b]:
            if f == i: continue
            o = bef(f, i)
            for c in labs:
                # after f the x_b wire carries Rx|Rt at f
                m.AddBoolOr([o.Not(), Rx[c, f].Not(), Rx[c, i]])
                m.AddBoolOr([o.Not(), Rt[c, f].Not(), Rx[c, i]])
        for f in bytgt[a]:
            if f == i: continue
            o = bef(f, i)
            for c in labs:
                m.AddBoolOr([o.Not(), Rx[c, f].Not(), Rt[c, i]])
                m.AddBoolOr([o.Not(), Rt[c, f].Not(), Rt[c, i]])
    if hint:
        for i, t in enumerate(hint): m.AddHint(T[i], t)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = time; s.parameters.num_workers = workers
    st = s.Solve(m)
    name = s.StatusName(st)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        order = sorted(range(E), key=lambda i: s.Value(T[i]))
        assert replay(B, edges, order, x0, t0) is None
        return 'LAMINAR', order
    return name, None


def sym_edges(pairs):
    out = []
    for a, b in pairs:
        out += [(a, b), (b, a)]
    return out


if __name__ == '__main__':
    # abstract free cliques: B = I (all distinct pairs free)
    for k in (2, 3, 4):
        B = [[int(i == j) for j in range(k)] for i in range(k)]
        E = sym_edges(list(itertools.combinations(range(k), 2)))
        r = brute(B, E, limit=10**6) if len(E) <= 8 else None
        print('free clique K%d bidirected: brute' % k, 'none' if r is None else r, ' sat:', sat(B, E, time=30)[0])
    # one-way (antisymmetric) clique corrections: only a->b for a<b
    for k in (3, 4, 5, 6):
        B = [[int(i == j) for j in range(k)] for i in range(k)]
        E = [(a, b) for a, b in itertools.combinations(range(k), 2)]
        print('free clique K%d one-way (tournament):' % k, sat(B, E, time=30)[0])
    # star: source 0 read by k targets that are pairwise required (B=1 between targets), bidirected
    for k in (2, 3, 5):
        n = k + 1
        B = [[1 if i == j or (i and j) else 0 for j in range(n)] for i in range(n)]
        E = sym_edges([(0, j) for j in range(1, n)])
        print('bidirected star, leaves pairwise required, k=%d:' % k, sat(B, E, time=30)[0])
        B2 = [[int(i == j) for j in range(n)] for i in range(n)]
        print('bidirected star, leaves pairwise free,     k=%d:' % k, sat(B2, E, time=30)[0])
    print('PASS_LAM')
