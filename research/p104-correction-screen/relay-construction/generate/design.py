"""Affine-halving relay network design for the P104 X-only correction (heuristic search).
Pools are linear frames V of F2^9; a pool's demand D(V) is a space of forms over sources in V.
Output: design.json with pools, children assignment and per-target split."""
import json, sys, random, collections, itertools
from fractions import Fraction as Q
D = json.load(open('geo.json')); n = D['n']; q = D['q']; nt = D['Hnormal']; G = D['G']
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 0; MODE = sys.argv[2] if len(sys.argv) > 2 else 'greedy'
random.seed(SEED)
GZ = list(range(1, 512)); random.Random(SEED).shuffle(GZ)
dot = lambda a, b: bin(a & b).count('1') & 1
def frame_of(basis):
    s = {0}
    for b in basis:
        if b not in s: s |= {x ^ b for x in s}
    m = 0
    for v in s: m |= 1 << v
    return m
def members(m): return [v for v in range(512) if m >> v & 1]
def dimf(m): return (m.bit_count()).bit_length() - 1
ustar = [u for u in range(1, 512) if all(dot(u, x) for x in q)]; assert len(ustar) == 1; US = ustar[0]
qmask = {s: 1 << q[s] for s in range(n)}
def srcs(m): return [s for s in range(n) if m >> q[s] & 1]
def perp_in(m, z):
    r = 0
    for v in members(m):
        if not dot(v, z): r |= 1 << v
    return r
# C = I - G/2 ; forms are dicts s->Fraction
C = [{j: (Q(1) if i == j else Q(0)) - Q(G[i][j], 2) for j in range(n) if (Q(1) if i == j else Q(0)) - Q(G[i][j], 2)} for i in range(n)]
def restrict(f, S): return {s: f[s] for s in S if f.get(s)}
def rref(vecs):
    """return list of basis dicts (reduced) of span"""
    basis = []  # (pivot, vec)
    for v in vecs:
        v = dict(v)
        for p, b in basis:
            if v.get(p):
                c = v[p]
                for k, x in b.items():
                    y = v.get(k, 0) - c*x
                    if y: v[k] = y
                    else: v.pop(k, None)
        if v:
            p = min(v); c = v[p]; v = {k: x/c for k, x in v.items()}
            nb = []
            for pp, b in basis:
                if b.get(p):
                    cc = b[p]; b = dict(b)
                    for k, x in v.items():
                        y = b.get(k, 0) - cc*x
                        if y: b[k] = y
                        else: b.pop(k, None)
                nb.append((pp, b))
            basis = nb + [(p, v)]
    return basis
def rank(vecs): return len(rref(vecs))
Hm = [perp_in((1 << 512) - 1, nt[t]) for t in range(n)]
FULL = (1 << 512) - 1
# ---------- top level: choose split per target ----------
def halves(V):
    """parallel-class splits of pool V: list of (Z1, Z2) (both hyperplanes of V, Z1∩Z2 ⊂ US-perp)"""
    V0 = perp_in(V, US); seen = set(); out = []
    for z in range(1, 512):
        Z = perp_in(V, z)
        if Z == V or Z == V0 or Z in seen: continue
        Z2 = perp_in(V, z ^ US); seen |= {Z, Z2}; out.append((Z, Z2))
    return out
demand = collections.defaultdict(list)       # pool -> list of forms (requests from parents)
requests = []                                  # (child, parent, form) edges: each is one lifetime move carrying a basis form
top = {}
splits_t = {t: halves(Hm[t]) for t in range(n)}
# choose target splits maximizing K-sharing via a matching-like greedy (pairs with proportional restriction first)
choice = {}
if MODE.startswith('ilp'):
    from ortools.sat.python import cp_model
    m = cp_model.CpModel(); x = {(t, i): m.NewBoolVar('') for t in range(n) for i in range(len(splits_t[t]))}
    for t in range(n): m.AddExactlyOne([x[t, i] for i in range(len(splits_t[t]))])
    halfuse = collections.defaultdict(list)
    for t in range(n):
        for i, (Z1, Z2) in enumerate(splits_t[t]):
            for Z in (Z1, Z2): halfuse[Z].append((t, i))
    obj = []
    for Z, us in halfuse.items():
        ts = sorted({t for t, _ in us})
        if len(ts) < 2: continue
        assert len(ts) == 2
        a, b = ts; S = srcs(Z); r = rank([restrict(C[a], S), restrict(C[b], S)])
        y = m.NewBoolVar(''); m.AddBoolOr([x[a, i] for tt, i in us if tt == a]).OnlyEnforceIf(y)
        for tt in (a, b): m.Add(sum(x[tt, i] for t2, i in us if t2 == tt) >= 1).OnlyEnforceIf(y)
        obj.append((3 - r) * y * 10 + y)   # prefer proportional shares
    m.Maximize(sum(obj)); sol = cp_model.CpSolver(); sol.parameters.max_time_in_seconds = 60; sol.parameters.num_workers = 2
    sol.Solve(m)
    for t in range(n):
        choice[t] = next(i for i in range(len(splits_t[t])) if sol.Value(x[t, i]))
else:
    for t in range(n): choice[t] = random.randrange(len(splits_t[t]))
for t in range(n):
    Z1, Z2 = splits_t[t][choice[t]]
    for Z in (Z1, Z2):
        f = restrict(C[t], srcs(Z))
        if f: demand[Z].append(f); requests.append((Z, ('H', t), f))
# ---------- shrink: replace a pool by the span of its sources when smaller ----------
def shrink(V):
    S = srcs(V); b = []
    for s in S:
        v = q[s]
        for x in b: v = min(v, v ^ x)
        if v: b.append(v)
    return frame_of(b)
# recursive levels
def process(level_pools, choose):
    nxt = collections.defaultdict(list); reqs = []; info = {}
    for V in level_pools:
        forms = demand[V]; S = srcs(V)
        if dimf(V) == 1:
            continue
        cands = halves(V)
        if 'glob' in MODE:
            V0 = perp_in(V, US); gc = []
            for z in GZ:
                Z = perp_in(V, z)
                if Z != V and Z != V0:
                    gc = [(Z, perp_in(V, z ^ US))]; break
            cands = gc or cands
        best = None
        for (Z1, Z2) in cands:
            cost = 0; parts = []
            for Z in (Z1, Z2):
                SZ = [s for s in S if Z >> q[s] & 1]
                if not SZ: continue
                Zs = shrink(Z)
                rf = rref([restrict(f, SZ) for f in forms]); r = len(rf)
                old = rank(demand[Zs]) if Zs in demand else 0
                new = rank(demand.get(Zs, []) + [b for _, b in rf])
                cost += r + (new - old) + (0.5 if Zs not in demand and Zs not in nxt else 0) + 0.01 * dimf(Zs)
                parts.append((Zs, [b for _, b in rf]))
            key = (cost + (random.random() * choose if choose else 0))
            if best is None or key < best[0]: best = (key, parts)
        for Zs, fs in best[1]:
            for f in fs:
                demand[Zs].append(f); reqs.append((Zs, V, f))
            nxt[Zs]
        info[V] = best[1]
    return reqs
# process pools in decreasing dimension order (children always smaller)
done = set(); allreq = list(requests)
while True:
    pend = [V for V in demand if V not in done and dimf(V) > 1]
    if not pend: break
    d = max(dimf(V) for V in pend)
    lvl = [V for V in pend if dimf(V) == d]
    lvl.sort(key=lambda V: -rank(demand[V]))
    allreq += process(lvl, 0.3 if MODE.endswith('rand') else 0)
    done |= set(lvl)
# accounting: per pool entrants (requests whose parent is V), exits (requests whose child is V)
ent = collections.Counter(); ex = collections.Counter()
for ch, par, f in allreq:
    ex[ch] += 1
    if par[0] != 'H' if isinstance(par, tuple) else True: ent[par] += 1
pools = set(demand)
births = {V: max(0, ex[V] - ent[V]) for V in pools}
L = sum(births.values()); visits = sum(max(ex[V], ent[V]) for V in pools)
lines = sum(1 for V in pools if dimf(V) == 1)
print('seed', SEED, MODE, 'pools', len(pools), 'L', L, 'visits', visits, 'requests', len(allreq), 'linepools', lines,
      'byDim', sorted(collections.Counter(dimf(V) for V in pools).items()))
json.dump(dict(L=L, visits=visits, choice=choice, mode=MODE, seed=SEED,
               req=[[str(ch), (par if isinstance(par, tuple) else str(par)), {str(k): str(v) for k, v in f.items()}] for ch, par, f in allreq]),
          open('design_%s_%d.json' % (MODE, SEED), 'w'))
