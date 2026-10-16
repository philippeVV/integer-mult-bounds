"""Phase-frame (H0-type) cost model reduced to F2 polar forms ("beta model").

A phase frame C_Q = H^{(x)m} diag(i^{Q(x)}) H^{(x)m}, Q a Z/4-valued quadratic form on F2^m.
Modulo free corrections (constant phases, translations X_l = 2[l.x] terms, Z-sign conjugations)
Q is determined by its polar form beta(Q) = B mod 2 in Sym_m(F2).  A kernel (aI+bX_v)^{+-1} adds
+-[v.x], i.e. adds v v^T to beta.  Hence the exact edge cost is
    c(A) = min #{v_j : A = sum v_j v_j^T} = rank_F2(A)        if A = 0 or diag(A) != 0
                                          = rank_F2(A) + 1    if A != 0 alternating,
and the endpoint condition is beta_out(rho w) - beta_in(w) = I (polar form of wt(x)).
Symmetric matrices are stored as bitmasks over the upper triangle (i<=j)."""
import json, os, itertools

S = os.path.dirname(os.path.abspath(__file__))
G = json.load(open(os.path.join(S, 'bh_graph.json')))   # copy of expN1/bh_graph.json (BH scenario A)
PHYS = [tuple(e) for e in G['edges']]            # 735 gate-gate segments (orientation irrelevant)
SESS = [tuple(s) for s in G['sessions']]          # 157 (source gate, sink gate)
NV = 182


def pairs(m):
    return [(i, j) for i in range(m) for j in range(i, m)]


def to_mat(mask, m):
    P = pairs(m); A = [[0] * m for _ in range(m)]
    for b, (i, j) in enumerate(P):
        if mask >> b & 1: A[i][j] = A[j][i] = 1
    return A


def rank2(A):
    rows = [int(''.join(map(str, r)), 2) for r in A]; r = 0
    for bit in reversed(range(len(A))):
        piv = next((k for k in range(r, len(rows)) if rows[k] >> bit & 1), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for k in range(len(rows)):
            if k != r and rows[k] >> bit & 1: rows[k] ^= rows[r]
        r += 1
    return r


def tables(m):
    P = pairs(m); nb = len(P)
    rk, cost = [], []
    for mask in range(1 << nb):
        A = to_mat(mask, m); r = rank2(A)
        alt = mask != 0 and all(A[i][i] == 0 for i in range(m))
        rk.append(r); cost.append(r + (1 if alt else 0))
    I = sum(1 << b for b, (i, j) in enumerate(P) if i == j)
    return nb, rk, cost, I


def min_vvT(A_mask, m):
    """brute-force check of c(A) for small m: min number of v v^T summing to A."""
    P = pairs(m)
    outer = []
    for v in range(1, 1 << m):
        msk = 0
        for b, (i, j) in enumerate(P):
            if (v >> i & 1) and (v >> j & 1): msk |= 1 << b
        outer.append(msk)
    if A_mask == 0: return 0
    for t in range(1, m + 3):
        for combo in itertools.combinations_with_replacement(outer, t):
            x = 0
            for o in combo: x ^= o
            if x == A_mask: return t
    return None


def s_of(beta, cost, I):
    s = sum(cost[beta[u] ^ beta[v]] for u, v in PHYS)
    s += sum(cost[beta[a] ^ beta[b] ^ I] for a, b in SESS)
    return s


if __name__ == '__main__':
    for m in (1, 2, 3):
        nb, rk, cost, I = tables(m)
        bad = [x for x in range(1 << nb) if min_vvT(x, m) != cost[x]]
        print(f"m={m}: |Sym|={1<<nb}, c-table verified by brute force: {not bad}")
