"""expXW: bootstrap (Boot) ceilings for NON-half-rate MDS local codes on PG(2,q), q <= 9.

Non-reversible networks: a gate is any linear map F^{in} -> F^{out}; zero sources and erasers are
dangling half-edges (cost 0, dimension-neutral).  On the internal graph the local code at v is the
graph of a map with k_v = in-degree, so a rate profile k is exactly an in-degree vector of an
orientation of Gamma (sum k = E).  For MDS local codes of any dims (repetition [d,1], parity [d,d-1],
...) the expXE bootstrap theorem gives dim C <= Boot(Gamma,k) = min_sigma sum_v (k_v - j_v)^+.
Every greedy value printed is an UPPER bound on Boot, hence on s = dim C / E.

Parts
  A  sanity: exact Boot by subset DP on PG(2,2); greedy agrees; reversal symmetry Boot(k)=Boot(d-k);
     Boot(k) <= FAS(any orientation with in-degrees k)  (checked on the constructive orientation)
  B  class-uniform profiles (t, d-t), t=1..d-1, all q
  C  mixed rep/parity profiles: fraction phi of points repetition (k=1) and of lines parity (k=q),
     rest near-half; first-N and random choice of the special vertices
  D  random orientations (k = in-degrees)
  E  hill-climb: maximise greedy-Boot over orientations by single-edge flips (adversarial for us)
Thresholds: s_thr(g=6, a) = 1/(1+5^{1-a});  a=1/9 -> 0.1930.
usage: python3 -I -B xw.py [restarts]   (stdlib only)
"""
import sys, random, math, json, itertools

def diffset(q):
    """Perfect difference set mod q^2+q+1 by backtracking (q <= 9)."""
    N = q * q + q + 1
    best = [0, 1]
    def rec(D, used):
        if len(D) == q + 1: return D
        for x in range(D[-1] + 1, N):
            new = set()
            ok = True
            for d in D:
                for y in ((x - d) % N, (d - x) % N):
                    if y in used or y in new: ok = False; break
                    new.add(y)
                if not ok: break
            if ok:
                r = rec(D + [x], used | new)
                if r: return r
        return None
    D = rec([0, 1], {1, N - 1})
    assert D and len(D) == q + 1
    return D

def pg(q):
    """Incidence graph of PG(2,q): vertices 0..N-1 points, N..2N-1 lines; edge list."""
    N = q * q + q + 1; edges = []
    if q > 9:   # use expXE's Singer construction (galois-built) instead of backtracking; limits-2 copies only SG16
        import os
        G = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'graphs_sg16.json')))[f'SG{q}']
        ends = [[] for _ in range(G['E'])]
        for v, es in enumerate(G['Pts'] + G['Lns']):
            for e in es: ends[e].append(v)
        edges = [tuple(x) for x in ends]
        adj = [[] for _ in range(2 * N)]
        for a, b in edges: adj[a].append(b); adj[b].append(a)
        assert all(len(x) == q + 1 for x in adj)
        return N, edges, adj
    D = diffset(q)
    for j in range(N):
        for d in D: edges.append(((j + d) % N, N + j))
    adj = [[] for _ in range(2 * N)]
    for a, b in edges: adj[a].append(b); adj[b].append(a)
    assert all(len(x) == q + 1 for x in adj)
    return N, edges, adj

def greedy(adj, k, rng):
    """Bucket greedy: always place an unplaced vertex of minimum residual (k_v - j_v)^+."""
    n = len(adj); j = [0] * n; placed = [False] * n
    K = max(k); buckets = [[] for _ in range(K + 1)]; where = [None] * n
    def put(w):
        d = max(k[w] - j[w], 0); b = buckets[d]; where[w] = (d, len(b)); b.append(w)
    def take(w):
        d, i = where[w]; b = buckets[d]; last = b[-1]; b[i] = last; where[last] = (d, i); b.pop()
    for w in range(n): put(w)
    c = 0; order = []
    for _ in range(n):
        d = next(d for d in range(K + 1) if buckets[d])
        b = buckets[d]; v = b[rng.randrange(len(b))]; take(v); placed[v] = True; c += d; order.append(v)
        for w in adj[v]:
            if not placed[w]:
                take(w); j[w] += 1; put(w)
    return c, order

def cost(adj, k, order):
    pos = [0] * len(adj)
    for i, v in enumerate(order): pos[v] = i
    return sum(max(k[v] - sum(1 for w in adj[v] if pos[w] < pos[v]), 0) for v in order)

def boot(adj, k, R, seed=1):
    rng = random.Random(seed); best = None
    for _ in range(R):
        c, o = greedy(adj, k, rng)
        if best is None or c < best[0]: best = (c, o)
    assert cost(adj, k, best[1]) == best[0]
    return best[0]

def exact_boot(adj, k):
    n = len(adj); INF = 10 ** 9; f = [INF] * (1 << n); f[0] = 0
    nb = [sum(1 << w for w in adj[v]) for v in range(n)]
    for S in range(1 << n):
        if f[S] == INF: continue
        for v in range(n):
            if not S >> v & 1:
                c = f[S] + max(k[v] - bin(S & nb[v]).count('1'), 0)
                T = S | 1 << v
                if c < f[T]: f[T] = c
    return f[(1 << n) - 1]

def orient(edges, n, k):
    """Return an orientation (list of (tail, head)) with in-degree k, or None (path-reversal method)."""
    if sum(k) != len(edges): return None
    arcs = [(a, b) for a, b in edges]          # a -> b
    indeg = [0] * n
    for a, b in arcs: indeg[b] += 1
    out = [[] for _ in range(n)]
    for i, (a, b) in enumerate(arcs): out[a].append(i)
    while True:
        under = [v for v in range(n) if indeg[v] < k[v]]
        if not under: break
        # BFS from an under-full w along directed arcs to an over-full u, then reverse the path
        w = under[0]; prev = {w: None}; Q = [w]; hit = None
        for x in Q:
            if indeg[x] > k[x]: hit = x; break
            for i in out[x]:
                y = arcs[i][1]
                if y not in prev: prev[y] = i; Q.append(y)
        if hit is None: return None
        x = hit
        while prev[x] is not None:
            i = prev[x]; a, b = arcs[i]
            out[a].remove(i); arcs[i] = (b, a); out[b].append(i); x = a
        indeg[w] += 1; indeg[hit] -= 1
    return arcs

def fas_upper(arcs, n, rng, R=50):
    """Upper bound on the minimum feedback arc set: best of R greedy (Eades-Lin-Smyth-like) orders."""
    best = None
    for _ in range(R):
        ind = [0] * n; outd = [0] * n
        for a, b in arcs: outd[a] += 1; ind[b] += 1
        alive = set(range(n)); outs = [[] for _ in range(n)]; ins = [[] for _ in range(n)]
        for a, b in arcs: outs[a].append(b); ins[b].append(a)
        left, right = [], []
        while alive:
            v = max(alive, key=lambda x: (outd[x] - ind[x], rng.random()))
            left.append(v); alive.discard(v)
            for w in outs[v]:
                if w in alive: ind[w] -= 1
            for w in ins[v]:
                if w in alive: outd[w] -= 1
        pos = {v: i for i, v in enumerate(left)}
        f = sum(1 for a, b in arcs if pos[a] > pos[b])
        best = f if best is None else min(best, f)
    return best

def thr(a, g=6): return 1 / (1 + (g - 1) ** (1 - a))
def a_of(beta, g=6):
    return -math.inf if beta <= 0 else 1 - math.log((1 - beta) / beta) / math.log(g - 1)

def main():
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    T = thr(1 / 9); out = {'threshold_a1/9': T}
    print(f'threshold s_thr(g=6,a=1/9) = {T:.4f}; a=0 -> {thr(0):.4f}; a=1/2 -> {thr(0.5):.4f}')
    # ---- A: sanity on PG(2,2) ----
    N, edges, adj = pg(2); E = len(edges)
    for t in (1, 2):
        k = [t] * N + [3 - t] * N
        ex = exact_boot(adj, k); gr = boot(adj, k, 400)
        print(f'A PG(2,2) t={t}: exact Boot {ex}, greedy {gr}, E {E}')
        assert ex == gr
    rng = random.Random(5)
    for trial in range(20):
        arcs = [(a, b) if rng.random() < .5 else (b, a) for a, b in edges]
        k = [0] * (2 * N)
        for a, b in arcs: k[b] += 1
        kr = [3 - x for x in k]
        e1, e2 = exact_boot(adj, k), exact_boot(adj, kr)
        assert e1 == e2, (e1, e2)
        assert e1 <= fas_upper(arcs, 2 * N, rng)
    print('A reversal symmetry Boot(k)=Boot(d-k) and Boot<=FAS: 20 random orientations of PG(2,2) PASS')
    res = {}
    for q in (2, 3, 4, 5, 7, 8, 9):
        N, edges, adj = pg(q); E = len(edges); d = q + 1; n = 2 * N
        row = {'E': E}
        # ---- B: class-uniform ----
        B = {}
        for t in range(1, d):
            k = [t] * N + [d - t] * N
            b = boot(adj, k, R); B[t] = b
        row['class_uniform'] = B
        bB = max(B.values()); tB = max(B, key=B.get)
        print(f'q={q} E={E} B class-uniform Boot by t: ' + ' '.join(f'{t}:{B[t]}' for t in B) +
              f'  max {bB} (t={tB}) s<={bB/E:.4f}')
        # ---- C: mixed rep/parity ----
        lo, hi = d // 2, d - d // 2           # near-half base: points lo, lines hi
        C = {}
        for phi in (0, .25, .5, .75, 1):
            m = round(phi * N); vals = []
            for mode in ('first', 'rand0', 'rand1', 'rand2'):
                if mode == 'first': P = list(range(m)); L = list(range(m))
                else:
                    r = random.Random(f'{q}|{phi}|{mode}')   # limits-2: deterministic (hash() of a str is salted per process)
                    P = r.sample(range(N), m); L = r.sample(range(N), m)
                k = [lo] * N + [hi] * N
                for p in P: k[p] = 1
                for l in L: k[N + l] = q
                if sum(k) != E: continue
                if orient(edges, n, k) is None: vals.append(None); continue
                vals.append(boot(adj, k, R))
            C[phi] = vals
        row['mixed_rep_parity'] = {str(p): v for p, v in C.items()}
        print('   C mixed (phi: first,rand x3): ' + '  '.join(f'{p}:{v}' for p, v in C.items()))
        # ---- D: random orientations ----
        rr = random.Random(q); Dv = []
        while len(Dv) < 12:
            k = [0] * n
            for a, b in edges:
                k[b if rr.random() < .5 else a] += 1
            if min(k) < 1 or max(k) > d - 1: continue   # pure sources/erasers are dangling, not gates
            Dv.append(boot(adj, k, R))
        row['random_orient'] = Dv
        print(f'   D random orientations: max {max(Dv)} mean {sum(Dv)/len(Dv):.1f}')
        # ---- E: hill-climb (maximise greedy Boot) ----
        if E <= 500:
            k = [lo] * N + [hi] * N
            arcs = orient(edges, n, k); cur = boot(adj, k, max(R // 3, 10), 11); hr = random.Random(99)
            steps = 300 if E < 200 else 150
            for it in range(steps):
                i = hr.randrange(E); a, b = arcs[i]
                k2 = list(k); k2[b] -= 1; k2[a] += 1
                if min(k2) < 1 or k2[a] > len(adj[a]) - 1: continue
                c2 = boot(adj, k2, max(R // 3, 10), 11)
                if c2 >= cur:
                    k, cur = k2, c2; arcs[i] = (b, a)
            final = boot(adj, k, 4 * R, 3)
            hist = {}
            for x in k: hist[x] = hist.get(x, 0) + 1
            row['hill_climb'] = dict(boot=final, k_hist=hist)
            print(f'   E hill-climb (greedy-Boot maximised over orientations): {final} s<={final/E:.4f} k-hist {sorted(hist.items())}')
        allv = list(B.values()) + [x for v in C.values() for x in v if x is not None] + Dv + \
               ([row['hill_climb']['boot']] if 'hill_climb' in row else [])
        mx = max(allv); row['max'] = mx; row['max_s'] = mx / E; row['a_ceiling'] = a_of(mx / E)
        row['beats_threshold'] = mx / E > T
        print(f'   => max Boot over all profiles {mx}, s<={mx/E:.4f}, a<={a_of(mx/E):.3f}, '
              f'{"ABOVE" if mx/E > T else "below"} threshold {T:.4f}')
        res[q] = row
    out['results'] = res
    json.dump(out, open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'xw.json'), 'w'), indent=1)
    if not any(r['beats_threshold'] for r in res.values()): print('NO PROFILE BEATS THRESHOLD on PG(2,q<=9)')
    print('PASS_EXPXW')

if __name__ == '__main__':
    main()
