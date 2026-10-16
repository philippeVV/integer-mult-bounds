#!/usr/bin/env python3
"""expQY stdlib checker: two-port GL(2,F_p) circuits that compute a permutation (monomial) matrix -- token routability.

  python3 qy_check.py frame  [trials]   Lemma 1 (token-frame reduction) on random F2 circuits
  python3 qy_check.py line   [trials]   Theorem 2 (nearest-neighbour line => routable) over F2 and F3, exact MITM
  python3 qy_check.py c2gap            the subset/LGV test (C2) is necessary but NOT sufficient for routing,
                                       and the C2-passing non-routable rho of the example is F2-unrealisable
  python3 qy_check.py hunt W n trials seed   C2-prefiltered exact MITM hunt for a counterexample (F2)
  python3 qy_check.py all              frame + line + c2gap + a small hunt

Conventions: a circuit is a list of (a, b, g); g = (p, q, r, s) means out_a = p*a + q*b, out_b = r*a + s*b (mod P).
The state v[pos] is the vector (tuple over inputs) carried by wire pos. rho maps input -> output wire.
Token routing: every gate independently passes or swaps its two tokens; S(C) = set of resulting rho.
"""
import itertools, random, sys

# ---------------------------------------------------------------- basic objects
def gl2(P):
    return [(p, q, r, s) for p in range(P) for q in range(P) for r in range(P) for s in range(P) if (p * s - q * r) % P]

def apply(v, a, b, g, P):
    p, q, r, s = g; va, vb = v[a], v[b]; v = list(v)
    v[a] = tuple((p * x + q * y) % P for x, y in zip(va, vb))
    v[b] = tuple((r * x + s * y) % P for x, y in zip(va, vb))
    return tuple(v)

def inv2(g, P):
    p, q, r, s = g; d = pow((p * s - q * r) % P, P - 2, P)
    return (s * d % P, -q * d % P, -r * d % P, p * d % P)

def ident(W):
    return tuple(tuple(int(i == j) for j in range(W)) for i in range(W))

def run(W, C, P):
    v = ident(W)
    for a, b, g in C: v = apply(v, a, b, g, P)
    return v

def monomial_perm(W, v):
    """rho (input -> wire) if the state is a monomial matrix, else None."""
    rho = [None] * W
    for pos, vec in enumerate(v):
        nz = [w for w in range(W) if vec[w]]
        if len(nz) != 1 or rho[nz[0]] is not None: return None
        rho[nz[0]] = pos
    return tuple(rho)

def routable_set(W, pairs):
    st = {tuple(range(W))}                         # st[pos] = token
    for a, b in pairs:
        nx = set(st)
        for s in st:
            l = list(s); l[a], l[b] = l[b], l[a]; nx.add(tuple(l))
        st = nx
    out = set()
    for s in st:
        r = [0] * W
        for p, t in enumerate(s): r[t] = p
        out.add(tuple(r))
    return out

def realisable_perms(W, pairs, P):
    """All rho such that some choice of GL(2,F_P) gates on `pairs` computes a monomial matrix D*P_rho (exact MITM)."""
    G = gl2(P); k = len(pairs) // 2
    F = {ident(W)}
    for a, b in pairs[:k]: F = {apply(v, a, b, g, P) for v in F for g in G}
    Bk = set()
    for rho in itertools.permutations(range(W)):
        for D in itertools.product(range(1, P), repeat=W):
            t = [None] * W
            for w in range(W): t[rho[w]] = tuple(D[w] if j == w else 0 for j in range(W))
            Bk.add((tuple(t), rho))
    out = set()
    # backward: states X with (remaining gates) X = target  <=>  X = g_k^-1 ... target
    cur = {}
    for t, rho in Bk: cur.setdefault(t, set()).add(rho)
    for a, b in reversed(pairs[k:]):
        nx = {}
        for t, rhos in cur.items():
            for g in G:
                u = apply(t, a, b, inv2(g, P), P)
                nx.setdefault(u, set()).update(rhos)
        cur = nx
    for t, rhos in cur.items():
        if t in F: out |= rhos
    return out

# ---------------------------------------------------------------- Lemma 1: token frame
PASS_TYPE = {(1, 0, 0, 1): 'id', (1, 1, 0, 1): 'a+=b', (1, 0, 1, 1): 'b+=a'}
TAU = (0, 1, 1, 0)

def mul(g, h, P=2):   # matrix of "h then g"
    p, q, r, s = g; A, B, Cc, D = h
    return ((p * A + q * Cc) % P, (p * B + q * D) % P, (r * A + s * Cc) % P, (r * B + s * D) % P)

def token_frame(W, C, x):
    """Given routing x (x_k = 1: swap), return (pi, Ctilde): tokens' positions pi (token -> wire) after the routing and the
    circuit on TOKEN labels whose k-th gate is pi_{k-1}^{-1} (tau^{x_k} g_k) pi_{k-1}."""
    pos = list(range(W)); tok = list(range(W))     # pos[token], tok[wire]
    Ct = []
    for (a, b, g), xk in zip(C, x):
        h = mul(TAU, g) if xk else g               # tau^{x_k} * g_k  (tau applied after g)
        ta, tb = tok[a], tok[b]
        Ct.append((ta, tb, h))                     # same 2x2 action, now on token labels (ta at "a" slot)
        if xk:
            tok[a], tok[b] = tb, ta; pos[ta], pos[tb] = b, a
    return tuple(pos), Ct

def cmd_frame(trials=300):
    G = gl2(2); bad = 0
    for t in range(trials):
        W = random.randint(2, 5); n = random.randint(1, 9)
        C = [(*random.sample(range(W), 2), random.choice(G)) for _ in range(n)]
        x = [random.randint(0, 1) for _ in range(n)]
        pi, Ct = token_frame(W, C, x)
        # (i) matrices: C = P_pi * Ctilde  (as states: wire pi(w) of C carries what token-slot w of Ctilde carries)
        v, vt = run(W, C, 2), run(W, Ct, 2)
        ok1 = all(v[pi[w]] == vt[w] for w in range(W))
        # (ii) routing sets: S(C) = pi o S(Ctilde)
        S = routable_set(W, [(a, b) for a, b, _ in C])
        St = routable_set(W, [(a, b) for a, b, _ in Ct])
        ok2 = S == {tuple(pi[s[w]] for w in range(W)) for s in St}
        # (iii) with x = pass/swap type of each gate, every Ctilde gate is the identity or a transvection
        xs = [0 if g in PASS_TYPE else 1 for _, _, g in C]
        ok3 = all(h in PASS_TYPE for _, _, h in token_frame(W, C, xs)[1])
        bad += not (ok1 and ok2 and ok3)
    print(f'[frame] Lemma 1 on {trials} random F2 circuits: failures {bad}')
    return bad == 0

# ---------------------------------------------------------------- Theorem 2: line architecture
def cmd_line(trials=60):
    bad = 0; tot = 0
    for P, Wmax, nmax in ((2, 5, 8), (3, 3, 5)):
        for t in range(trials):
            W = random.randint(3, Wmax); n = random.randint(2, nmax)
            pairs = [(i, i + 1) for i in (random.randrange(W - 1) for _ in range(n))]
            R = realisable_perms(W, pairs, P); S = routable_set(W, pairs); tot += 1
            if not R <= S: bad += 1; print('  LINE VIOLATION', P, W, pairs, R - S)
            assert S <= R      # routings are realisable (pure pass/swap gates)
    print(f'[line] Theorem 2 (adjacent pairs, F2 and F3, monomial targets, exact MITM): {tot} circuits, violations {bad}')
    return bad == 0

# ---------------------------------------------------------------- C2 is necessary, not sufficient
def c2_candidates(W, pairs):
    R = routable_set(W, pairs); im = {}
    for A in range(1, (1 << W) - 1):
        st = {A}
        for a, b in pairs:
            nx = set(st)
            for s in st:
                if ((s >> a) ^ (s >> b)) & 1: nx.add(s ^ (1 << a) ^ (1 << b))
            st = nx
        im[A] = st
    return [rho for rho in itertools.permutations(range(W)) if rho not in R and
            all(sum(1 << rho[w] for w in range(W) if A >> w & 1) in S for A, S in im.items())]

def cmd_c2gap():
    W = 4; pairs = [(0, 2), (1, 2), (1, 3), (0, 3), (0, 2), (0, 1)]
    c = c2_candidates(W, pairs); R = realisable_perms(W, pairs, 2)
    ok = bool(c) and not (set(c) & R) and R == routable_set(W, pairs)
    print(f'[c2gap] pairs {pairs}: C2-passing non-routable rho = {c}; F2-realisable among them: {sorted(set(c) & R)}; '
          f'realisable set == routable set: {R == routable_set(W, pairs)}')
    return ok

# ---------------------------------------------------------------- hunt
def cmd_hunt(W, n, trials, seed):
    random.seed(seed); allp = list(itertools.combinations(range(W), 2)); nc = 0; found = []
    for t in range(trials):
        pairs = [random.choice(allp) for _ in range(n)]
        c = c2_candidates(W, pairs)
        if not c: continue
        nc += len(c); R = realisable_perms(W, pairs, 2)
        found += [(pairs, rho) for rho in c if rho in R]
    print(f'[hunt] W={W} n={n} trials={trials}: C2-candidates {nc}, realisable non-routable {len(found)}', found[:3])
    return not found

if __name__ == '__main__':
    random.seed(1)
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if cmd == 'frame': ok = cmd_frame(int(sys.argv[2]) if len(sys.argv) > 2 else 300)
    elif cmd == 'line': ok = cmd_line(int(sys.argv[2]) if len(sys.argv) > 2 else 60)
    elif cmd == 'c2gap': ok = cmd_c2gap()
    elif cmd == 'hunt': ok = cmd_hunt(*map(int, sys.argv[2:6]))
    else: ok = all([cmd_frame(), cmd_line(), cmd_c2gap(), cmd_hunt(4, 6, 400, 1), cmd_hunt(5, 8, 200, 2)])
    print('ALL OK' if ok else 'FAILED'); print('PASS_EXPQY_CHECK' if ok else 'FAIL_EXPQY_CHECK'); sys.exit(0 if ok else 1)
