"""Free graphs for two non-plane-generic label alphabets on the 112 selected ports:
 (A) paired planes (PR1825 matching; E_i = E_j = <q_i,q_j>) in the rank-1 term model;
 (B) F8-line ports: F2^9 = F8^3 with x.y = Tr(sum x_i y_i) (self-dual basis), port label L_a = F8*q_a (3 unpaid).
Reports edges, clique number bound, alpha (CP-SAT, certified optimal), and the floors that follow."""
import json, random, sys
from common import *
from ortools.sat.python import cp_model

P = load_ports(); n = len(P)


def alpha(edges, n, tl=120):
    m = cp_model.CpModel(); x = [m.NewBoolVar('x%d' % i) for i in range(n)]
    for a, b in edges:
        m.Add(x[a] + x[b] <= 1)
    m.Maximize(sum(x)); s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl
    s.parameters.num_workers = 2
    st = s.Solve(m)
    return int(s.ObjectiveValue()), int(s.BestObjectiveBound()), s.StatusName(st)


res = {}
# (A) paired planes
cert = json.load(open(MAIN + '/research/paired-plane-port-matching-seed6c563a/CERTIFICATE.json'))
plane = {}
for i, j in cert['matching']:
    plane[i] = plane[j] = (P[i], P[j])
edges = [(a, b) for a in range(n) for b in range(a + 1, n)
         if all(dot(x, y) == 0 for x in plane[a] for y in plane[b])]
assert all(not (set(plane[a]) == set(plane[b])) for a, b in edges)
al = alpha(edges, n)
res['paired'] = dict(free_edges=len(edges), alpha=al)
print('paired planes: free edges', len(edges), 'alpha', al, '-> any decoder C >= alpha+1; PSD (cost>=2/term) C >= 2*alpha')

# (B) F8 structure
MOD = 0b1011  # t^3+t+1
def mul(a, b):
    r = 0
    for i in range(3):
        if b >> i & 1: r ^= a << i
    for i in (4, 3):
        if r >> i & 1: r ^= MOD << (i - 3)
    return r
def tr(x): return (x ^ mul(x, x) ^ mul(mul(x, x), mul(x, x))) & 7
assert all(tr(x) in (0, 1) for x in range(8))
sd = None
import itertools
for b in itertools.combinations(range(1, 8), 3):
    if all(tr(mul(b[i], b[j])) == (i == j) for i in range(3) for j in range(3)):
        if rank_f2(list(b)) == 3: sd = b; break
assert sd, 'no self-dual basis'
print('self-dual basis of F8/F2:', sd)


def structure(perm):
    """phi(x)_i = sum over coords k with perm[k] in block i of x_k * sd[perm[k]%3]."""
    def phi(x):
        v = [0, 0, 0]
        for k in range(9):
            if x >> k & 1:
                v[perm[k] // 3] ^= sd[perm[k] % 3]
        return v
    return phi


def ortho_mix(x, vs):
    for v in vs:
        if dot(x, v): x ^= v
    return x


def f8_graph(perm, vs):
    phi = structure(perm)
    F = [phi(ortho_mix(q, vs)) for q in P]
    def h(u, v):
        s = 0
        for a, b in zip(u, v): s ^= mul(a, b)
        return s
    # sanity: trace form reproduces the dot product (on the transformed labels)
    T = [ortho_mix(q, vs) for q in P]
    assert all(tr(h(F[a], F[b])) == dot(T[a], T[b]) for a in range(n) for b in range(n))
    assert all(dot(T[a], T[b]) == dot(P[a], P[b]) for a in range(n) for b in range(n))
    free = [(a, b) for a in range(n) for b in range(a + 1, n) if h(F[a], F[b]) == 0]
    # same F8 line?
    def line(u):
        return frozenset(tuple(mul(l, c) for c in u) for l in range(1, 8))
    L = [line(u) for u in F]
    shared = sum(1 for a in range(n) for b in range(a + 1, n) if L[a] == L[b])
    return free, shared, len(set(L))

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
rows = []
trials = int(sys.argv[2]) if len(sys.argv) > 2 else 12
for t in range(trials):
    perm = list(range(9)); random.shuffle(perm)
    vs = [] if t == 0 else [random.choice([v for v in range(1, 512) if bin(v).count('1') % 2 == 0]) for _ in range(6)]
    free, shared, nl = f8_graph(perm if t else list(range(9)), vs)
    al = alpha(free, n, tl=60)
    rows.append(dict(trial=t, free_edges=len(free), distinct_lines=nl, shared_pairs=shared, alpha=al))
    print('F8 trial', t, 'free edges', len(free), 'distinct F8 lines', nl, 'port pairs sharing a line', shared, 'alpha', al)
res['f8'] = rows
json.dump(res, open('b_graphs.json', 'w'), indent=1)
