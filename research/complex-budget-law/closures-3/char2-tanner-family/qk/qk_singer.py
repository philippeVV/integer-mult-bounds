"""expQK: Singer-block dim W for the Hermitian Tanner code on PG(2,q) - (Singer matching), any q = 2^e.
Port of expQB/qb.py (same construction, same gauge) with log-table arithmetic for F_{q^3}, so q = 64
(F_{2^18}) is feasible.  Only edges at i = 0 (plus a sample) are built explicitly; the cocycle theta is
read from sampled wrap edges and asserted consistent.
usage: qk_singer.py e j m0   (e=2 j=1 -> expect 9; e=4 j=2 -> expect 772)"""
import sys, time, random, collections, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # gf.py (expNC)
from gf import GF, IRR

PRIMPOLY = {9: (1 << 9) | (1 << 4) | 1, 15: (1 << 15) | 2 | 1, 6: (1 << 6) | 0b11, 12: (1 << 12) | (1 << 6) | (1 << 4) | 2 | 1, 18: (1 << 18) | (1 << 7) | 1}

class LF:
    """F_{2^m} via log/exp tables (no n x n table)."""
    def __init__(s, m, poly):
        n = 1 << m; s.n = n; s.ORD = n - 1
        EXP = np.zeros(2 * s.ORD, dtype=np.int64); LOG = np.full(n, -1, dtype=np.int64); v = 1
        for k in range(s.ORD):
            EXP[k] = EXP[k + s.ORD] = v; LOG[v] = k; v <<= 1
            if v & n: v ^= poly
        assert v == 1 and (LOG[1:] >= 0).all(), "poly not primitive"
        s.EXP, s.LOG = EXP, LOG
    def mul(s, a, b):
        a = np.asarray(a, dtype=np.int64); b = np.asarray(b, dtype=np.int64)
        la = s.LOG[a]; lb = s.LOG[b]
        r = s.EXP[(la + lb) % s.ORD]
        return np.where((a == 0) | (b == 0), 0, r)
    def inv(s, a): return int(s.EXP[(-s.LOG[a]) % s.ORD])
    def pw(s, a, k):
        a = int(a)
        if a == 0: return 1 if k == 0 else 0
        return int(s.EXP[(int(s.LOG[a]) * k) % s.ORD])
    def rank(s, M):
        M = np.array(M, dtype=np.int64); r = 0; nr, nc = M.shape
        for c in range(nc):
            if r >= nr: break
            nz = np.nonzero(M[r:, c])[0]
            if len(nz) == 0: continue
            i = r + nz[0]
            if i != r: M[[r, i]] = M[[i, r]]
            M[r] = s.mul(s.inv(int(M[r, c])), M[r])
            col = M[:, c].copy(); col[r] = 0; rows = np.nonzero(col)[0]
            if len(rows): M[rows] ^= s.mul(col[rows][:, None], M[r][None, :])
            r += 1
        return r

def model(e, j, m0=0, verbose=True, nsample=400):
    q = 2 ** e; k = q // 2; N = q * q + q + 1
    F = GF(2, e); B = LF(3 * e, PRIMPOLY[3 * e])
    FL = np.full(q, -1, dtype=np.int64); FE = np.zeros(2 * (q - 1), dtype=np.int64)
    # F_q log tables
    gq = next(g for g in range(2, q) if len({F.pw(g, t) for t in range(q - 1)}) == q - 1) if q > 2 else 1
    v = 1
    for t in range(q - 1): FE[t] = FE[t + q - 1] = v; FL[v] = t; v = int(F.mul[v, gq])
    def fpw(a, n):
        if a == 0: return 1 if n == 0 else 0
        return int(FE[(FL[a] * n) % (q - 1)])
    irr = IRR[(2, e)]
    def ev(th, cs):
        r = 0; p = 1
        for c in cs:
            if c: r ^= p
            p = int(B.mul(p, th))
        return r
    th = next(t for t in range(2, B.n) if ev(t, list(irr)) == 0)
    phi = np.zeros(q, dtype=np.int64)
    for a in range(q):
        vv = 0; p = 1
        for i in range(e):
            if a >> i & 1: vv ^= p
            p = int(B.mul(p, th))
        phi[a] = vv
    for a in range(q):
        for b in range(q): assert phi[F.mul[a, b]] == B.mul(phi[a], phi[b])
    g = 2; g2 = int(B.mul(g, g)); g3 = int(B.mul(g2, g))
    # solve phi[c0] + phi[c1] g + phi[c2] g^2 = g3 : linear over F_2 -> brute over q^3 is 2.6e5 for q=64, ok
    PG1 = B.mul(phi, g); PG2 = B.mul(phi, g2)
    sol = []
    for c2 in range(q):
        t = (PG2[c2] ^ g3)
        for c1 in range(q):
            t2 = t ^ PG1[c1]
            hit = np.nonzero(phi == t2)[0]
            for c0 in hit: sol.append((int(c0), c1, c2))
    assert len(sol) == 1; k0, k1, k2 = sol[0]
    A = np.array([[0, 0, k0], [1, 0, k1], [0, 1, k2]], dtype=np.int16)
    def mv(Mx, v):
        out = [0, 0, 0]
        for r in range(3):
            for c in range(3): out[r] ^= int(F.mul[Mx[r, c], v[c]])
        return tuple(out)
    def minv(Mx):
        aug = np.hstack([Mx, np.eye(3, dtype=np.int16)]); R, piv = F.rref(aug); assert piv == [0, 1, 2]
        return R[:, 3:].copy()
    AT = A.T.copy(); AmT = minv(AT)
    def dot(x, y):
        r = 0
        for a, b in zip(x, y): r ^= int(F.mul[a, b])
        return r
    def norm(t):
        for x in t:
            if x: iv = int(F.inv[x]); return tuple(int(F.mul[iv, y]) for y in t)
    p0 = (1, 0, 0)
    pencil = [(0, 0, 1)] + [(0, 1, a) for a in range(q)]
    l0 = pencil[m0]
    P = [p0]; Lr = [l0]
    for i in range(1, N): P.append(mv(A, P[-1])); Lr.append(mv(AmT, Lr[-1]))
    assert len({norm(x) for x in P}) == N and len({norm(x) for x in Lr}) == N
    assert norm(mv(A, P[-1])) == norm(p0)
    Dall = [d for d in range(N) if dot(p0, Lr[d]) == 0]; assert len(Dall) == q + 1 and Dall[0] == 0
    Dp = Dall[1:]; Ds = set(Dall)
    rng = random.Random(1)
    for i in rng.sample(range(N), min(N, 40)):
        for jj in range(N): assert (dot(P[i], Lr[jj]) == 0) == (((jj - i) % N) in Ds)
    def coords(u, v, y):
        # solve lam*u + mu*v = y over F_q
        for a, b in [(0, 1), (0, 2), (1, 2)]:
            det = int(F.add[F.mul[u[a], v[b]], F.mul[u[b], v[a]]])
            if det: break
        di = int(F.inv[det])
        lam = int(F.mul[di, F.add[F.mul[y[a], v[b]], F.mul[y[b], v[a]]]])
        mu = int(F.mul[di, F.add[F.mul[u[a], y[b]], F.mul[u[b], y[a]]]])
        for c in range(3): assert int(F.add[F.mul[lam, u[c]], F.mul[mu, v[c]]]) == y[c]
        return lam, mu
    def mono(lam, mu): return [int(F.mul[fpw(lam, k - 1 - r), fpw(mu, r)]) for r in range(k)]
    fr = lambda x: F.frob(x, j)
    def lperp(l):
        out = []
        for x in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
            if dot(x, l) == 0: out.append(x)
        # solve generically: two independent solutions of l.x = 0
        for a in range(q):
            for b in range(q):
                for c in range(q):
                    x = (a, b, c)
                    if any(x) and dot(x, l) == 0:
                        out.append(x)
                        if len(out) >= q + 3: return out
        return out
    u0, v0 = (0, 1, 0), (0, 0, 1)
    cand = lperp(l0); u0p = cand[0]; v0p = next(x for x in cand if norm(x) != norm(u0p))
    # bases: Pu[i] = A^{-T i} u0, Lu[j] = A^j u0p  -> build iteratively (cheap)
    Pu = [u0]; Pv = [v0]; Lu = [u0p]; Lv = [v0p]
    for i in range(1, N):
        Pu.append(mv(AmT, Pu[-1])); Pv.append(mv(AmT, Pv[-1])); Lu.append(mv(A, Lu[-1])); Lv.append(mv(A, Lv[-1]))
    def edge(i, d):
        jj = (i + d) % N
        lam, mu = coords(Pu[i], Pv[i], Lr[jj]); gg = mono(lam, mu)
        lam, mu = coords(Lu[jj], Lv[jj], P[i]); hh = mono(fr(lam), fr(mu))
        return gg, hh
    E0 = {d: edge(0, d) for d in Dp}
    G0 = {d: E0[d][0] for d in Dp}; H0 = {d: E0[d][1] for d in Dp}
    def ratio(x, y):
        nzi = next(r for r in range(k) if y[r]); s = int(F.mul[x[nzi], F.inv[y[nzi]]])
        assert all(int(F.mul[s, y[r]]) == x[r] for r in range(k)); return s
    theta = None; nwrap = 0
    samp = rng.sample(range(N), min(N, nsample)) + list(range(N - max(Dp), N))
    for i in samp:
        for d in Dp:
            gg, hh = edge(i, d)
            s = ratio(gg, G0[d]); t = ratio(hh, H0[d]); rt = int(F.mul[t, F.inv[s]])
            if i + d < N: assert rt == 1, (i, d, s, t)
            else:
                nwrap += 1
                if theta is None: theta = rt
                assert rt == theta
    if theta is None: theta = 1
    lt = int(B.LOG[phi[theta]]); assert lt % N == 0; rho = int(B.EXP[lt // N])
    Gm = np.array([[phi[x] for x in G0[d]] for d in Dp], dtype=np.int64)
    Hm = np.array([[int(B.mul(B.pw(rho, d), phi[x])) for x in H0[d]] for d in Dp], dtype=np.int64)
    zeta = int(B.EXP[q - 1])
    fac = [pp for pp in range(2, N + 1) if N % pp == 0 and all(pp % r for r in range(2, int(pp ** .5) + 1))]
    assert B.pw(zeta, N) == 1 and all(B.pw(zeta, N // pp) != 1 for pp in fac)
    if verbose: print(f"[model] q={q} j={j} m0={m0} N={N} theta={theta} wrap-edges checked={nwrap}", flush=True)
    return dict(B=B, q=q, k=k, N=N, Dp=Dp, G=Gm, H=Hm, zeta=zeta)

def spectrum(Mo):
    B = Mo['B']; N = Mo['N']; Dp = Mo['Dp']; G = Mo['G']; H = Mo['H']; z = Mo['zeta']; q = Mo['q']
    spec = []
    for chi in range(N):
        lam = np.array([B.pw(z, (chi * d) % N) for d in Dp], dtype=np.int64)
        M = np.hstack([G, B.mul(lam[:, None], H)])
        spec.append(q - B.rank(M))
    return spec

if __name__ == "__main__":
    e, j, m0 = map(int, sys.argv[1:4])
    t0 = time.time(); Mo = model(e, j, m0)
    print(f"model built {time.time()-t0:.1f}s", flush=True)
    sp = spectrum(Mo); S = sum(sp); q = Mo['q']; N = Mo['N']; E = N * q
    import math
    sig = S / E; mk = (E - S) / S
    a = 1 - math.log(mk) / math.log(5)
    print(f"q={q} j={j} m0={m0}: S=dimW={S} spectrum={dict(sorted(collections.Counter(sp).items()))} E={E} "
          f"sigma=S/E={sig:.4f} (E-S)/S={mk:.4f} a(b=5)={a:+.4f}  ({time.time()-t0:.1f}s)", flush=True)
    r = math.isqrt(q); closed = q * (q * q - r ** 3 + 1) // 4      # S(q) = q(q^2 - q^1.5 + 1)/4 (expQK)
    print(f"closed form S(q) = {closed}")
    print('PASS_EXPQK_SINGER' if S == closed else 'FAIL_EXPQK_SINGER')
