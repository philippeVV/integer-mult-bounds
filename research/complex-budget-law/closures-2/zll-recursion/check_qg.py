#!/usr/bin/env python3
"""expQG checker (stdlib only): Lemma W / W' (cycle form), tensor invariance, and the Lemma-W pre-filter.

Run:  python3 -I -B check_qg.py [--repo PATH]   (default: the repository containing this file)
Prints PASS_EXPQG at the end.

Part 1  Lemma W  (s >= m*|Fix rho|) and Lemma W' (s >= m*cyc(rho)) on random circuits, scored by the
        UPSTREAM contract code scripts/finite_bit_contract.py (physical_graph + inspect_candidate), plus a
        greedy adversarial frame minimiser that tries to push s below the bound.
Part 2  Tensor invariance of the tilt root: Psi_{S(x)S'} = Psi_S * Psi_S', root of a tensor lies in
        [min a, max a] and equals the secant-weighted average; the expQC 'mediant' is exact only for
        single-rank levels.
Part 3  Pre-filter on the live candidate schedules (inputs/schedule784.json from expNI, inputs/causal_*.json and
        psa_*.json from expNC, converted from pickle) and on ZLL compute-copy-uncompute compiles
        (inputs/bh_compile.json from expNN, cross-check).
"""
import json, math, os, random, sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
INPUTS = os.path.join(HERE, 'inputs')
REPO = os.environ.get('IMB_REPO') or os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
if '--repo' in sys.argv:
    REPO = sys.argv[sys.argv.index('--repo') + 1]
sys.path.insert(0, os.path.join(REPO, 'scripts'))
import finite_bit_contract as fbc          # upstream definitions of rho, physical segments, s
from audit_joint_frames import eye, rank, sub


def cycles(rho):
    seen, out = set(), []
    for w in range(len(rho)):
        if w in seen:
            continue
        c, x = [], w
        while x not in seen:
            seen.add(x); c.append(x); x = rho[x]
        out.append(c)
    return out


# ---------------------------------------------------------------- Part 1
def rand_circuit(rng, W):
    gates = []
    for _ in range(rng.randint(1, 5)):
        k = rng.randint(2, min(3, W))
        roles = rng.sample(range(W), k)
        xors = []
        for _ in range(rng.randint(1, 4)):
            t, s_ = rng.sample(roles, 2)
            xors.append([t, s_])
        gates.append(dict(roles=roles, xors=xors))
    return gates


def rand_mat(rng, m, lo=-1, hi=1):
    return [[rng.randint(lo, hi) for _ in range(m)] for _ in range(m)]


def score(W, m, gates, rho, src, gfr):
    I = eye(m)
    snk = [None] * W
    for w in range(W):
        snk[rho[w]] = [[Fr(x) for x in row] for row in
                       (tuple(a + b for a, b in zip(r1, r2)) for r1, r2 in zip(src[w], I))]
    cand = dict(W=W, m=m, gates=[dict(g, frame=f) for g, f in zip(gates, gfr)], rho=list(rho),
                source_frames=src, sink_frames=snk)
    res = fbc.inspect_candidate(cand, require_deficit=False)      # upstream scorer
    return res


def strand_sums(W, gates, ranks):
    edges = fbc.physical_graph(W, gates)
    tot = [0] * W
    for (u, v, role), r in zip(edges, ranks):
        tot[role] += r
    return tot


def switch_perms(W, gates):
    """All pi realisable by choosing one port permutation per gate (straight/cross 'switch' choices)."""
    from itertools import permutations, product
    opts = [list(permutations(range(len(g['roles'])))) for g in gates]
    for choice in product(*opts):
        pos = list(range(W))                 # pos[w] = role currently holding input w's strand content
        where = {w: w for w in range(W)}     # role -> input index held
        for g, p in zip(gates, choice):
            roles = g['roles']
            held = [where[r] for r in roles]
            for i, r in enumerate(roles):
                where[roles[p[i]]] = held[i]
        pi = [None] * W
        for r, w in where.items():
            pi[w] = r
        yield tuple(pi)


def switch_cycle_bound(W, gates, rho):
    rinv = [0] * W
    for w, r in enumerate(rho):
        rinv[r] = w
    best = 0
    for pi in switch_perms(W, gates):
        best = max(best, len(cycles([rinv[pi[w]] for w in range(W)])))
    return best


def part1(trials=1500, seed=1):
    print('== Part 1: Lemma W and W\' against the upstream contract scorer ==')
    rng = random.Random(seed)
    nperm = nfix = ntight = nbetween = nroute = 0
    worst = None
    for _ in range(trials):
        W = rng.randint(2, 5); m = rng.choice([2, 3])
        gates = rand_circuit(rng, W)
        try:
            rho = fbc.scalar_permutation(W, gates)
        except ValueError:
            continue
        nperm += 1
        cyc = cycles(rho)
        src = [rand_mat(rng, m) for _ in range(W)]
        gfr = [rand_mat(rng, m) for _ in gates]
        res = score(W, m, gates, rho, src, gfr)
        tot = strand_sums(W, gates, res['edge_ranks'])
        for c in cyc:                       # W' per cycle; W per fixed point (|c| = 1)
            assert sum(tot[w] for w in c) >= m, ('cycle strand sum < m', c, tot)
        nfix += sum(1 for c in cyc if len(c) == 1)
        assert res['s'] >= m * len(cyc)
        sc = switch_cycle_bound(W, gates, rho)          # pi = id gives cyc(rho^-1) = cyc(rho)
        assert sc >= len(cyc) and res['s'] >= m * sc
        nbetween += (len(cyc) < sc < W); nroute += (sc == W)
        # greedy adversary: coordinate descent over a small frame alphabet, trying to minimise s
        alph = [[[0] * m for _ in range(m)]] + [[[int(i == j) * t for j in range(m)] for i in range(m)]
                                                for t in (1, -1)] + \
               [[[int(i == j == k) for j in range(m)] for i in range(m)] for k in range(m)]
        best = res['s']
        for _ in range(2):
            for kind in ('g', 's'):
                arr = gfr if kind == 'g' else src
                for i in range(len(arr)):
                    keep = arr[i]
                    for A in alph:
                        arr[i] = A
                        sv = score(W, m, gates, rho, src, gfr)['s']
                        if sv < best:
                            best, keep = sv, A
                    arr[i] = keep
        assert best >= m * sc
        ratio = Fr(best, m * len(cyc))
        if best == m * len(cyc):
            ntight += 1
        if worst is None or ratio < worst[0]:
            worst = (ratio, W, m, len(cyc), best)
    print(f'  {nperm} permutation circuits (W<=5, m in 2,3), {nfix} fixed roles: every own strand and every '
          f'rho-cycle strand bundle has rank sum >= m; adversary attains s = m*cyc(rho) on {ntight} circuits')
    print(f'  switch-cycle bound s >= m*max_pi cyc(rho^-1 pi): holds everywhere (also vs adversary); '
          f'{nroute} circuits switch-routable (bound = Wm), {nbetween} strictly between cyc(rho) and W')
    # controls showing which hypotheses are load-bearing
    # (i) terminal segments must be charged: gate-free identity role, zero source, identity sink
    W, gates = 2, [dict(roles=[0, 1], xors=[[0, 1], [0, 1]])]           # x0 ^= x1 twice = identity
    rho = fbc.scalar_permutation(W, gates)
    edges = fbc.physical_graph(W, gates)
    term = [e for e in edges if e[0] < W or e[1] >= W + len(gates)]
    assert rho == (0, 1) and len(term) == 2 * W
    res = score(W, 2, gates, rho, [[[0, 0], [0, 0]]] * 2, [[[0, 0], [0, 0]]])
    interior = sum(r for (u, v, _), r in zip(edges, res['edge_ranks'])
                   if not (u < W or v >= W + len(gates)))
    print(f'  control (i): physical_graph charges all {len(term)} terminal segments; zero-frame identity '
          f'circuit s = {res["s"]} = Wm, but its interior-only charge is {interior} (bound needs terminals)')
    # (ii) cycle form needs char 0 (or char not dividing cycle length): over F2, 2*I = 0
    print('  control (ii): the cycle form telescopes to len(c)*I_m; rank m over Q, but 0 over F_p when p | len(c).'
          ' Fixed-point form (len 1) holds over every field.')
    return nperm


# ---------------------------------------------------------------- Part 2
def Psi(hist, m, W, a):
    return sum(n * (r / m) ** (1 - a) for r, n in hist.items()) / W


def root(hist, m, W):
    lo, hi = 0.0, 1.0
    if Psi(hist, m, W, 0.0) >= 1:
        return 0.0
    for _ in range(200):
        mid = (lo + hi) / 2
        (lo, hi) = (mid, hi) if Psi(hist, m, W, mid) < 1 else (lo, mid)
    return lo


def tensor(A, B):
    (h1, m1, W1), (h2, m2, W2) = A, B
    h = {}
    for r1, n1 in h1.items():
        for r2, n2 in h2.items():
            h[r1 * r2] = h.get(r1 * r2, 0) + n1 * n2
    return h, m1 * m2, W1 * W2


def dsum(A, B):
    (h1, m1, W1), (h2, m2, W2) = A, B
    assert m1 == m2
    h = dict(h1)
    for r, n in h2.items():
        h[r] = h.get(r, 0) + n
    return h, m1, W1 + W2


def rand_seed(rng):
    m = rng.randint(4, 30)
    W = rng.randint(3, 40)
    h = {}
    for _ in range(rng.randint(1, 4)):
        r = rng.randint(1, m - 1)
        h[r] = h.get(r, 0) + rng.randint(1, 3 * W)
    # rescale so that 0 < deficit: need sum r n < W m
    while sum(r * n for r, n in h.items()) >= W * m:
        W += 1
    return h, m, W


def part2(trials=2000, seed=2):
    print('\n== Part 2: tensor invariance of the tilt root ==')
    prof = {1: 349604, 2: 60192, 3: 62585, 4: 20448, 5: 11808, 7: 7, 8: 1152, 48: 4866}
    Wsel, a_star = 14336, 0.001675462676808
    msel = (sum(n * r ** (1 - a_star) for r, n in prof.items()) / Wsel) ** (1 / (1 - a_star))
    S1 = (prof, msel, Wsel)
    a1 = root(*S1); a2 = root(*tensor(S1, S1))
    print(f'  selected profile (m = {msel:.4f}): a(S) = {a1:.12f}, a(S(x)S) = {a2:.12f}')
    assert abs(a1 - a2) < 1e-12
    rng = random.Random(seed)
    worst_med = 0.0; n_ok = 0
    for _ in range(trials):
        A, B = rand_seed(rng), rand_seed(rng)
        aA, aB = root(*A), root(*B)
        T = tensor(A, B); aT = root(*T)
        for a in (0.0, 0.3, 0.77):
            assert abs(Psi(*T, a) - Psi(*A, a) * Psi(*B, a)) < 1e-9 * max(1, Psi(*T, a))
        assert min(aA, aB) - 1e-12 <= aT <= max(aA, aB) + 1e-12
        # exact secant-weight representation: a_T = (sA aA + sB aB)/(sA+sB), s_X = ln Psi_X(a_T)/(a_T - a_X)
        if abs(aT - aA) > 1e-9 and abs(aT - aB) > 1e-9:
            sA = math.log(Psi(*A, aT)) / (aT - aA); sB = math.log(Psi(*B, aT)) / (aT - aB)
            assert sA > 0 and sB > 0 and abs((sA * aA + sB * aB) / (sA + sB) - aT) < 1e-9
        assert abs(root(*tensor(A, A)) - aA) < 1e-9
        if A[1] == B[1]:
            aD = root(*dsum(A, B)); assert min(aA, aB) - 1e-12 <= aD <= max(aA, aB) + 1e-12
        # expQC mediant: sum ln(1/(1-delta_k)) / sum ln m_k  with delta = 1 - s/(Wm)
        med = sum(math.log(X[2] * X[1] / sum(r * n for r, n in X[0].items())) for X in (A, B)) / \
            (math.log(A[1]) + math.log(B[1]))
        worst_med = max(worst_med, abs(med - aT) / max(aT, 1e-12))
        n_ok += 1
    # single-rank levels: mediant (in the uniform metric with L = ln(m/r)) is exact
    A = ({3: 7}, 8, 3); B = ({2: 4}, 5, 2)
    aA, aB, aT = root(*A), root(*B), root(*tensor(A, B))
    LA, LB = math.log(8 / 3), math.log(5 / 2)
    med_exact = (aA * LA + aB * LB) / (LA + LB)
    assert abs(med_exact - aT) < 1e-12
    print(f'  {n_ok} random seed pairs: Psi multiplicative, a(A(x)B) in [min,max], equals the secant-weighted '
          f'average, a(A(x)A) = a(A); direct sums also in [min,max]')
    print(f'  expQC mediant sum ln(Wm/s)/sum ln m: worst relative error {worst_med:.2f} (NOT exact in general);'
          f' single-rank levels: (a_A L_A + a_B L_B)/(L_A+L_B) exact ({med_exact:.6f} = {aT:.6f})')


# ---------------------------------------------------------------- Part 3
def strands_from_schedule(pos, edges, R):
    """Orient: R edges backward, others forward in pos. Return (eulerian, arcs, back) ."""
    Rs = set(R)
    arcs = []
    for i, (x, y) in enumerate(edges):
        u, w = (x, y) if pos[x] < pos[y] else (y, x)
        arcs.append((w, u) if i in Rs else (u, w))
    V = len(pos)
    ind = [0] * V; outd = [0] * V
    for u, w in arcs:
        outd[u] += 1; ind[w] += 1
    imbalance = sum(abs(ind[v] - outd[v]) for v in range(V)) // 2
    return imbalance, arcs, Rs


def rho_from_pairing(arcs, back, succ):
    """register = start back arc e (input segment -> head(e)); follow succ through forward arcs to the
    next back arc f (its tail produces f -> output segment).  Returns end-map e -> f (cycle type of rho)."""
    end = {}
    for e in back:
        a = succ[e]; steps = 0
        while a not in back:
            a = succ[a]; steps += 1
            assert steps <= len(arcs)
        end[e] = a
    assert sorted(end.values()) == sorted(back)
    return end


def naive_pairing(arcs):
    ins, outs = {}, {}
    for i, (u, w) in enumerate(arcs):
        outs.setdefault(u, []).append(i); ins.setdefault(w, []).append(i)
    succ = {}
    for v in ins:
        for a, b in zip(sorted(ins[v]), sorted(outs[v])):
            succ[a] = b
    return succ


def euler_pairing(arcs):
    """Hierholzer: one closed trail per weakly connected component; succ along the circuit."""
    outs = {}
    for i, (u, w) in enumerate(arcs):
        outs.setdefault(u, []).append(i)
    ptr = {u: list(reversed(l)) for u, l in outs.items()}
    used = [False] * len(arcs); succ = {}; ncirc = 0
    for start in range(len(arcs)):
        if used[start]:
            continue
        ncirc += 1
        # iterative Hierholzer on arcs
        stack = [start]; used[start] = True; circuit = []
        while stack:
            a = stack[-1]; v = arcs[a][1]
            while ptr.get(v) and used[ptr[v][-1]]:
                ptr[v].pop()
            if ptr.get(v):
                b = ptr[v].pop(); used[b] = True; stack.append(b)
            else:
                circuit.append(stack.pop())
        circuit.reverse()
        for x, y in zip(circuit, circuit[1:] + circuit[:1]):
            assert arcs[x][1] == arcs[y][0]
            succ[x] = y
    return succ, ncirc


def cyc_of_map(mp):
    seen = set(); fix = 0; c = 0
    for e in mp:
        if e in seen:
            continue
        c += 1; x = e; L = 0
        while x not in seen:
            seen.add(x); x = mp[x]; L += 1
        fix += L == 1
    return fix, c


def part3():
    print('\n== Part 3: Lemma W / W\' pre-filter on live candidates ==')
    rows = []
    cands = [('expNI schedule784 (PG(2,16)-M, x^4, Singer d0=59)', os.path.join(INPUTS, 'schedule784.json'), 772)]
    for f, d in [('causal_q16_88_2_1_m0', 772), ('causal_q9_55_1_4_m-1', 166),
                 ('causal_q9_55_1_1_m-1', 166), ('psa_q9_55_1_1_m-1', 166),
                 ('causal_q8_54_1_1_m-1', 114)]:
        cands.append((f'expNC {f}', os.path.join(INPUTS, f + '.json'), d))
    for name, path, dimW in cands:
        if not os.path.exists(path):
            raise SystemExit(f'missing input {path}')
        d = json.load(open(path))
        if 'order' in d:
            pos = [0] * len(d['order'])
            for i, v in enumerate(d['order']):
                pos[v] = i
            edges = [tuple(a) for a in d['arcs']]; R = d['R']
        else:
            pos, R = d['pos'], d['R']
            N = len(pos) // 2                       # expNC causal.py: ends = (point, N + line)
            edges = [(pp, N + l) for (pp, l) in d['edges']]
        imb, arcs, back = strands_from_schedule(pos, edges, R)
        if imb:
            print(f'  {name}: |R| = {len(R)}, orientation not Eulerian (imbalance {imb}): registers are not '
                  f'conserved, needs C/J lanes -> not an all-role permutation circuit; Lemma W n/a')
            rows.append((name, len(R), dimW, None)); continue
        fn, cn = cyc_of_map(rho_from_pairing(arcs, back, naive_pairing(arcs)))
        succ, ncomp = euler_pairing(arcs)
        fe, ce = cyc_of_map(rho_from_pairing(arcs, back, succ))
        segs = len(arcs) + len(back)
        perm = 'permutation circuit' if len(R) == dimW else f'|R| > dim W = {dimW}: NOT a permutation (needs {len(R) - dimW} C+J)'
        print(f'  {name}: W = |R| = {len(R)} registers ({4 * len(R)} bit roles), {segs} segments/bit-lane; {perm}')
        print(f'      naive port pairing: |Fix| = {fn}, cyc = {cn} -> delta <= {1 - fn / len(R):.4f} (W), '
              f'{1 - cn / len(R):.4f} (W\')')
        print(f'      Euler-circuit pairing ({ncomp} component): |Fix| = {fe}, cyc = {ce} -> delta <= '
              f'{1 - fe / len(R):.4f} / {1 - ce / len(R):.4f}  => Lemma W cannot bite (pairing is free)')
        assert fe == 0 and ce == ncomp
        rows.append((name, len(R), dimW, (fn, cn, fe, ce)))
    # ZLL compute-copy-uncompute compiles: W = 3k + 2E roles; rho = identity except k S<->Y swaps
    print('  ZLL compute-copy-uncompute compiles (registers S,Y,p per session; sigma,tau per edge):')
    bhj = os.path.join(INPUTS, 'bh_compile.json')
    for name, k, E, Fstar in [('BH PG(2,9) (expNN)', 157, 735, 'F*=W (expNN, Lemma C)'),
                              ('PG(2,16)-M x^4 (hypothetical exact)', 772, 4368, 'F*=W if K*nu>=k (Lemma C)')]:
        W = 3 * k + 2 * E; fix = W - 2 * k; cyc = W - k
        print(f'    {name}: W = {W}, |Fix| = {fix}, cyc = {cyc}: delta <= 2k/W = {2 * k / W:.4f} (W), '
              f'k/W = {k / W:.4f} (W\'); binding: {Fstar}')
        if 'expNN' in name and os.path.exists(bhj):
            assert json.load(open(bhj))['W'] == W and json.load(open(bhj))['F_star'] == W
    print('    general: delta <= k/(3k+2E) <= 1/(3+2E/k)  (halves expQC\'s 2/(3+2m/k))')
    return rows


if __name__ == '__main__':
    part1()
    part2()
    part3()
    print('PASS_EXPQG')
