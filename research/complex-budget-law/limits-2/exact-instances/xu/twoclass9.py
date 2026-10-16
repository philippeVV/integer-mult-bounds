"""expXU copy of expXN/twoclass.py on the exact lane-free q = 9 schedule (only the import changed).
Two-class theorem check (all m, all X, Y at once).

Two-class terminal frames: beta_k = X for k in S, Y for k not in S (X, Y arbitrary in Sym_m(F2), any m).
Order sessions S | S^c.  Z = A0 (x) I + A1 (x) N with A0 = diag(T_SS, T_S'S'), A1 = offdiag(T_SS', T_S'S), N = I + X + Y.
Over F2-bar put N in Jordan form: det Z = prod_mu det(A0 + mu A1)^{mult}, mu over eigenvalues of N (mu in F_{2^d}, d <= m).
det(A0 + mu A1) = det T_SS det T_S'S' det(I - mu^2 M_S),  M_S = T_SS^{-1} T_SS' T_S'S'^{-1} T_S'S.
Hence: if, for ONE char-2 gate draw G, T_SS and T_S'S' are invertible and charpoly(M_S) has no root in F_{2^d}
for d = 1..D, then for EVERY m <= D and EVERY X, Y the GRS-F2 nullity is 0, i.e. s >= Wm for all interior frames.
Root test: deg gcd(chi_M(x), x^{2^d} - x) over GF(2^16) (detects roots in F_{2^d} inside F2-bar for any d).
usage: twoclass.py npartitions D
"""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import screen9 as S
g = S.grs2
EXP, LOG, Q = g.EXP, g.LOG, g.Q


def mul(a, b): return 0 if a == 0 or b == 0 else int(EXP[(int(LOG[a]) + int(LOG[b])) % (Q - 1)])
def inv(a): return int(EXP[(Q - 1 - int(LOG[a])) % (Q - 1)])


def ginv_mat(A):
    n = A.shape[0]; M = np.concatenate([A.copy(), np.eye(n, dtype=np.int64)], axis=1)
    for c in range(n):
        nz = np.nonzero(M[c:, c])[0]
        if len(nz) == 0: return None
        p = c + nz[0]; M[[c, p]] = M[[p, c]]
        M[c] = g.gmul(M[c], inv(int(M[c, c])))
        rows = np.nonzero(M[:, c])[0]; rows = rows[rows != c]
        if len(rows): M[rows] ^= g.gmul(M[rows, c][:, None], M[c][None, :])
    return M[:, n:]


def gdet(A):
    M = A.copy(); n = M.shape[0]; d = 1
    for c in range(n):
        nz = np.nonzero(M[c:, c])[0]
        if len(nz) == 0: return 0
        p = c + nz[0]; M[[c, p]] = M[[p, c]]
        piv = int(M[c, c]); d = mul(d, piv)
        M[c] = g.gmul(M[c], inv(piv))
        rows = np.nonzero(M[c + 1:, c])[0] + c + 1
        if len(rows): M[rows] ^= g.gmul(M[rows, c][:, None], M[c][None, :])
    return d


# ---- polynomials over GF(2^16): lists, low degree first ----
def ptrim(p):
    while p and p[-1] == 0: p.pop()
    return p
def padd(a, b):
    r = [0] * max(len(a), len(b))
    for i, x in enumerate(a): r[i] ^= x
    for i, x in enumerate(b): r[i] ^= x
    return ptrim(r)
def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1) if a and b else []
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y: r[i + j] ^= mul(x, y)
    return ptrim(r)
def pmod(a, f):
    a = a[:]; df = len(f) - 1; il = inv(f[-1])
    while len(ptrim(a)) - 1 >= df:
        c = mul(a[-1], il); sh = len(a) - 1 - df
        for i, y in enumerate(f):
            if y: a[sh + i] ^= mul(c, y)
        ptrim(a)
    return a
def pgcd(a, b):
    a, b = ptrim(a[:]), ptrim(b[:])
    while b: a, b = b, pmod(a, b)
    return a


def interp(xs, ys):
    """Lagrange interpolation over GF(2^16)."""
    n = len(xs); res = []
    for i in range(n):
        num = [1]; den = 1
        for j in range(n):
            if j != i: num = pmul(num, [xs[j], 1]); den = mul(den, xs[i] ^ xs[j])
        res = padd(res, [mul(c, mul(ys[i], inv(den))) for c in num])
    return res


def charpoly(M):
    n = M.shape[0]; xs = list(range(1, n + 2)); I = np.eye(n, dtype=np.int64)
    ys = [gdet(M ^ (I * x)) for x in xs]
    return interp(xs, ys)


def subfield_root_degrees(f, D):
    """list of d <= D with deg gcd(f, x^{2^d} - x) > 0."""
    out = []; xp = [0, 1]
    for d in range(1, D + 1):
        xp = pmod(pmul(xp, xp), f)                  # x^{2^d} mod f
        h = pgcd(f, padd(xp, [0, 1]))
        if len(h) > 1: out.append((d, len(h) - 1))
    return out


def check(T, Sset, D):
    W = T.shape[0]; Sc = [k for k in range(W) if k not in set(Sset)]
    A = T[np.ix_(Sset, Sset)]; Bm = T[np.ix_(Sset, Sc)]; C = T[np.ix_(Sc, Sset)]; Dm = T[np.ix_(Sc, Sc)]
    Ai, Di = ginv_mat(A), ginv_mat(Dm)
    if Ai is None or Di is None: return 'principal minor singular'
    M = g.gmatmul(g.gmatmul(Ai, Bm), g.gmatmul(Di, C))
    f = charpoly(M); v = 0
    while f[v] == 0: v += 1                          # eigenvalue 0 of M_S is irrelevant (lambda = 1/mu^2 != 0)
    return subfield_root_degrees(f[v:], D) if len(f) - v > 1 else []


if __name__ == '__main__':
    npart, D = int(sys.argv[1]), int(sys.argv[2])
    D = min(D, 15)                                   # entries live in GF(2^16): certify d <= 15
    W = S.W; rng = np.random.default_rng(4242)
    parts = {'src_parity': [k for k in range(W) if S.SRC[k] % 2],
             'src_point': [k for k in range(W) if S.SRC[k] < S.c.N],
             'half_B_order': list(range(W // 2))}
    for nm, key in (('snk_parity', lambda k: S.SNK[k] % 2), ('src_lt_snk_mod3', lambda k: (S.SRC[k] + S.SNK[k]) % 3 == 0)):
        parts[nm] = [k for k in range(W) if key(k)]
    for i in range(npart):
        parts[f'random{i}'] = sorted(rng.choice(W, size=int(rng.integers(20, W - 20)), replace=False).tolist())
    Ts = [g.transfer(np.random.default_rng(s)) for s in S.GSEEDS]
    res = {}
    for name, Sset in parts.items():
        r = [check(T, Sset, D) for T in Ts[:1]]
        if r[0]:                                    # roots found on draw 1: re-test on the other draws
            r += [check(T, Sset, D) for T in Ts[1:]]
        res[name] = {'|S|': len(Sset), 'roots_in_F2^d (d, multiplicity) per G draw': r}
        print(f"{name:14s} |S|={len(Sset):3d}  subfield roots of chi(M_S), d<= {D}: {r}", flush=True)
    json.dump(res, open(os.path.join(HERE, 'twoclass.json'), 'w'), indent=1)
