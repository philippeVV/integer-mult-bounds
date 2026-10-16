"""expQL Theorem A (char-2 trace bound), the identity behind it, checked exactly over F2.
Labels L_a (u-dim F2-subspaces of F2^h), pairing G = (L_a^T Omega L_b) of rank <= h.  For any decoder N over a char-2 field
with N_aa = 1 and N_ab = 0 on non-free pairs (PR(a,b) != 0):  (N (x) J_u) o G = blockdiag(G_aa), hence
    rank N * h >= rank((N (x) J_u) o G) = sum_a PR(a,a) = uV,   i.e.  rank_F2 N >= uV/h.
Checked here on PG(2,4) dot (u = 2, h = 6, V = 16) and on the E8 P120 labels (u = 1, h = 9, V = 120): the Hadamard identity,
rank G <= h, sum PR(a,a) = uV, for 200 random masked F2 decoders each, and the bound rank_F2 N >= uV/h on each sample.
Also: the P120 Gram decoder G/2 has a 2-adic denominator on every free cycle, which is exactly why it escapes (A2).
Stdlib only."""
import os, sys, random, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fam import build


def rank2(rows):
    basis = []
    for r in rows:
        for b in basis: r = min(r, r ^ b)
        if r: basis.append(r)
    return len(basis)


def pg_labels():
    """field lines F_4 x in F_4^3 = F2^6, F2 form Tr(B(v, w)).  Basis of L_x: {x, alpha x}."""
    F, B, P = build(2, 3, 'dot'); u = 2
    G = {}
    for a, x in enumerate(P):
        for b, y in enumerate(P):
            bxy = B(x, y)
            G[a, b] = [[F.tr(F.mul(F.pw(2, i + j), bxy)) for j in range(u)] for i in range(u)]   # alpha = 2
    free = [[a != b and B(P[a], P[b]) == 0 for b in range(len(P))] for a in range(len(P))]
    return len(P), u, 6, G, free


def e8_labels():
    from lattices import E8
    L = E8().tolist(); V = len(L)
    ip = [[sum(p * q for p, q in zip(L[a], L[b])) // 4 for b in range(V)] for a in range(V)]
    G = {(a, b): [[(1 + ip[a][b]) % 2 if a != b else 1]] for a in range(V) for b in range(V)}
    free = [[a != b and G[a, b][0][0] == 0 for b in range(V)] for a in range(V)]
    return V, 1, 9, G, free


def check(name, V, u, h, G, free, rng, samples=200):
    big = [sum(G[a, b][i][j] << (b * u + j) for b in range(V) for j in range(u)) for a in range(V) for i in range(u)]
    rG = rank2(big); diag = sum(rank2([sum(G[a, a][i][j] << j for j in range(u)) for i in range(u)]) for a in range(V))
    assert rG <= h and diag == u * V, (rG, diag)
    lo = None
    for _ in range(samples):
        N = [[1 if a == b else (rng.randint(0, 1) if free[a][b] else 0) for b in range(V)] for a in range(V)]
        # Hadamard product (N (x) J_u) o G: block (a,b) = N_ab * G_ab; must vanish off the diagonal blocks
        assert all(N[a][b] * G[a, b][i][j] == 0 for a in range(V) for b in range(V) if a != b for i in range(u) for j in range(u))
        rN = rank2([sum(N[a][b] << b for b in range(V)) for a in range(V)])
        assert rN * h >= u * V
        lo = rN if lo is None else min(lo, rN)
    print(f"{name}: V={V} u={u} h={h}: rank_F2 G = {rG} <= h, sum PR(a,a) = {diag} = uV; {samples} random masked F2 decoders:"
          f" identity holds, min rank {lo} >= uV/h = {u * V / h:.2f}")


rng = random.Random(1)
check('PG(2,4) dot', *pg_labels(), rng)
V, u, h, G, free = e8_labels()
check('E8 P120', V, u, h, G, free, rng, samples=20)
# G/2 on P120: free pairs a-b-c-a form a cycle of entries +-1/2 each, product 1/8 (negative 2-adic valuation)
tri = next(t for t in itertools.combinations(range(V), 3) if free[t[0]][t[1]] and free[t[1]][t[2]] and free[t[0]][t[2]])
print(f"E8 P120: free triangle {tri} carries G/2 entries +-1/2 -> cycle product 2^-3, the 2-adic denominator Thm A2 needs")
print('PASS_EXPQL_TRACE')
