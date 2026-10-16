"""Independent brute-force check (no CP-SAT, no rectangle lemma): enumerate every small graded class graph,
keep those satisfying the FULL entropy cut for every node set L, and test the locality theorem's claims.

usage: brute_check.py T p kmax mq md [nplant]
  T levels, p classes per level, k <= kmax extra climbs (labelled E-arcs inside levels),
  demand multiplicity <= mq and D multiplicity <= md per class pair.
Potential: union-find with offsets in Z^k (E_j arc u->v: c(v) = c(u) + e_j; demand: c(b) = c(a)).
Inconsistent E/demand data (an E-cycle, or two demands forcing different offsets) is discarded: such frames
cannot exist when I, E_1..E_k are Q-independent.
Claims checked on every surviving graph:
  C0  each D's endpoints lie in one (E + demand) component, so P_D = c(b) - c(a) is determined;
  C1  X(K,K') = 0 for every pair of E-components (K at level t, K' at level t+1);
  C2  all demands K -> K' share one class pair (hub);
  C3  supp P_D within E's of comp(a) u comp(b) (levels t, t+1), entries in {0, +-1};
  C4  each E_j has P_D[j] != 0 for <= 1 D (with multiplicity) per adjacent level pair (so <= 2 overall);
  C5  matrix plantings: E_j = s_j * (coordinate projection), m = 3; delta = sum_D (m - rank D) - R <= R.
"""
import itertools, sys, random
import numpy as np

T, p, kmax, mq, mdd = map(int, sys.argv[1:6])
nplant = int(sys.argv[6]) if len(sys.argv) > 6 else 4
n = T * p
node = lambda t, i: t * p + i
cells = [(t, a, b) for t in range(T - 1) for a in range(p) for b in range(p)]
nc = len(cells)
Eslots = [(t, u, v) for t in range(T) for u in range(p) for v in range(p) if u != v]
masks = np.arange(1, 1 << n)
inL = ((masks[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)  # (NL, n)
# cell enters L iff head in L and tail not in L
enterC = np.array([[inL[li, node(t + 1, b)] * (1 - inL[li, node(t, a)]) for (t, a, b) in cells]
                   for li in range(len(masks))], dtype=np.int64)  # (NL, nc)

def enterE(E):
    return np.array([sum(inL[li, node(t, v)] * (1 - inL[li, node(t, u)]) for (t, u, v) in E)
                     for li in range(len(masks))], dtype=np.int64)

def potentials(E, qsupp):
    """union-find with Z^k offsets; returns (root, off) or None if inconsistent"""
    k = len(E)
    par = list(range(n)); off = [np.zeros(k, dtype=np.int64) for _ in range(n)]
    def find(x):
        if par[x] == x: return x, np.zeros(k, dtype=np.int64)
        r, o = find(par[x]); return r, o + off[x]
    def union(x, y, delta):  # c(y) - c(x) = delta
        rx, ox = find(x); ry, oy = find(y)
        if rx == ry:
            return np.array_equal(oy - ox, delta)
        par[ry] = rx; off[ry] = ox + delta - oy; return True
    for j, (t, u, v) in enumerate(E):
        e = np.zeros(k, dtype=np.int64); e[j] = 1
        if not union(node(t, u), node(t, v), e): return None
    for ci in qsupp:
        (t, a, b) = cells[ci]
        if not union(node(t, a), node(t + 1, b), np.zeros(k, dtype=np.int64)): return None
    return [find(x) for x in range(n)]

def ecomp(E):
    par = list(range(n))
    def f(x):
        while par[x] != x: x = par[x]
        return x
    for (t, u, v) in E: par[f(node(t, u))] = f(node(t, v))
    return [f(x) for x in range(n)]

def rank(M):
    return np.linalg.matrix_rank(M.astype(float))

def multisets(supp_choices, mult, total=None):
    for vals in itertools.product(*[range(0 if not s else 1, mult + 1) if s else [0] for s in supp_choices]):
        yield vals

stats = dict(graphs=0, Dslots=0, unmatched=0, maxuse_pair=0, maxuse_total=0, max_delta_over_R=-9.0, fails=0)
rng = random.Random(1)
m = 3
for k in range(0, kmax + 1):
    for E in itertools.combinations(Eslots, k):
        # labelled E's: combinations suffice (labels are symmetric)
        if potentials(list(E), []) is None: continue  # E-cycle
        eE = enterE(E); comp = ecomp(E)
        # demand supports: subsets of cells, consistent
        for qs in range(1 << nc):
            qsupp = [i for i in range(nc) if qs >> i & 1]
            if not qsupp: continue
            pot = potentials(list(E), qsupp)
            if pot is None: continue
            for qm in itertools.product(range(1, mq + 1), repeat=len(qsupp)):
                qv = np.zeros(nc, dtype=np.int64); qv[qsupp] = qm
                qin = enterC @ qv  # demands entering each L
                per_t = [sum(qv[i] for i in range(nc) if cells[i][0] == t) for t in range(T - 1)]
                # D vectors with per-level totals equal (level cut equality is implied by cuts; enforce for speed)
                per_level_opts = []
                for t in range(T - 1):
                    idx = [i for i in range(nc) if cells[i][0] == t]
                    opts = [vals for vals in itertools.product(range(mdd + 1), repeat=len(idx)) if sum(vals) == per_t[t]]
                    per_level_opts.append((idx, opts))
                if any(not o for _, o in per_level_opts): continue
                for combo in itertools.product(*[o for _, o in per_level_opts]):
                    dv = np.zeros(nc, dtype=np.int64)
                    for (idx, _), vals in zip(per_level_opts, combo): dv[idx] = vals
                    if np.any(qin > enterC @ dv + eE): continue  # full entropy cut, every L
                    stats['graphs'] += 1
                    ok = True
                    # C0 / C3
                    Pd = {}
                    for i in np.nonzero(dv)[0]:
                        (t, a, b) = cells[i]
                        ra, oa = pot[node(t, a)]; rb, ob = pot[node(t + 1, b)]
                        if ra != rb: ok = False; print("C0 fail", E, qv, dv); break
                        P = ob - oa; Pd[i] = P
                        allowed = {j for j, (tt, u, v) in enumerate(E)
                                   if comp[node(tt, u)] in (comp[node(t, a)], comp[node(t + 1, b)])}
                        if any(P[j] != 0 and j not in allowed for j in range(k)) or np.any(np.abs(P) > 1):
                            ok = False; print("C3 fail", E, qv, dv, P)
                        stats['Dslots'] += int(dv[i])
                        if np.any(P): stats['unmatched'] += int(dv[i])
                    # C1 / C2
                    for t in range(T - 1):
                        Ks = {comp[node(t, a)] for a in range(p)}; Kps = {comp[node(t + 1, b)] for b in range(p)}
                        for K in Ks:
                            for Kp in Kps:
                                ids = [i for i in range(nc) if cells[i][0] == t and comp[node(t, cells[i][1])] == K
                                       and comp[node(t + 1, cells[i][2])] == Kp]
                                if sum(qv[i] - dv[i] for i in ids) != 0: ok = False; print("C1 fail", E, qv, dv)
                                if len([i for i in ids if qv[i]]) > 1: ok = False; print("C2 fail", E, qv, dv)
                    # C4
                    for j in range(k):
                        tot = 0
                        for t in range(T - 1):
                            u = sum(int(dv[i]) for i in Pd if cells[i][0] == t and Pd[i][j] != 0)
                            if u > 1: ok = False; print("C4 fail", E, qv, dv)
                            tot += u
                        stats['maxuse_pair'] = max(stats['maxuse_pair'], max([sum(int(dv[i]) for i in Pd if cells[i][0] == t and Pd[i][j] != 0) for t in range(T - 1)] or [0]))
                        stats['maxuse_total'] = max(stats['maxuse_total'], tot)
                    # C5 matrix plantings
                    if k and any(np.any(P) for P in Pd.values()):
                        for _ in range(nplant):
                            Em = []
                            for j in range(k):
                                supp = rng.sample(range(m), rng.randint(1, m))
                                Pj = np.zeros((m, m), dtype=np.int64)
                                for s_ in supp: Pj[s_, s_] = 1
                                Em.append(rng.choice([1, -1]) * Pj)
                            R = sum(rank(x) for x in Em)
                            dsum = sum(int(dv[i]) * (m - rank(np.eye(m, dtype=np.int64) + sum(int(Pd[i][j]) * Em[j] for j in range(k))))
                                       for i in Pd)
                            dr = (dsum - R) / R
                            stats['max_delta_over_R'] = max(stats['max_delta_over_R'], dr)
                            if dr > 1 + 1e-9: ok = False; print("C5 fail", E, dv, dr)
                    if not ok: stats['fails'] += 1
    print(f"T={T} p={p} k<={k} mq={mq} md={mdd}: {stats}", flush=True)
print('PASS_EXPRH_BRUTE' if stats['fails'] == 0 and stats['maxuse_total'] <= 1 and stats['max_delta_over_R'] <= 0 else 'FAIL_EXPRH_BRUTE')   # C0-C4 hold, single use, no deficit
