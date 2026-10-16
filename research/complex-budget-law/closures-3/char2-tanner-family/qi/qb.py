"""Singer-symmetric model of the PG(2,q) - (Singer perfect matching) twisted Tanner code, q = 2^e.
Coordinates over F_q (gf.py, same as expNC pg.py); Singer cycle A = multiplication by a primitive
element of F_{q^3} in the basis 1, g, g^2.  p_i = A^i p0, l_j = A^{-T j} l0 (i, j in [0,N)),
p0 = (1,0,0), l0 = m0-th line through p0; matching p_i <-> l_i.
All codes use the SAME global representatives p_i, l_j and the transported bases, so the code is
gauge-equivalent to expNC's (any global choice of representatives/bases is a gauge).
Wrap edges carry a cocycle theta; the F_{q^3}-gauge a_i -> rho^i a_i, b_j -> rho^j b_j (rho^N = theta)
plus row scalings turns the system into an exactly circulant one with rows (G_d, rho^d H_d)."""
import sys, numpy as np
import os; HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # gf.py (expNC)
from gf import GF, IRR

PRIM = {2: (6, 0b1000011), 3: (9, 0b1000010001), 4: (12, (1 << 12) | (1 << 6) | (1 << 4) | 2 | 1)}

class BF:
    """F_{2^m}, elements = bitmasks (poly basis, alpha = x)."""
    def __init__(s, m, poly):
        s.m, s.poly = m, poly; n = 1 << m; s.n = n; s.ORD = n - 1
        EXP = np.zeros(2 * s.ORD, dtype=np.int64); LOG = np.full(n, -1, dtype=np.int64); v = 1
        for k in range(s.ORD):
            EXP[k] = EXP[k + s.ORD] = v; LOG[v] = k; v <<= 1
            if v & n: v ^= poly
        assert v == 1 and len(set(EXP[:s.ORD].tolist())) == s.ORD
        s.EXP, s.LOG = EXP, LOG
        L = LOG[1:].astype(np.int32)
        MUL = np.zeros((n, n), dtype=np.int16)
        MUL[1:, 1:] = EXP[(L[:, None] + L[None, :]) % s.ORD].astype(np.int16)
        s.MUL = MUL
        s.INV = np.zeros(n, dtype=np.int16); s.INV[1:] = EXP[(-LOG[1:]) % s.ORD]
    def pw(s, a, k):
        a = int(a)
        if a == 0: return 1 if k == 0 else 0
        return int(s.EXP[(int(s.LOG[a]) * k) % s.ORD])
    def mul(s, a, b): return s.MUL[a, b]
    def rref(s, M):
        M = np.array(M, dtype=np.int16).copy(); r = 0; piv = []; nr, nc = M.shape; MUL = s.MUL
        for c in range(nc):
            if r >= nr: break
            nz = np.nonzero(M[r:, c])[0]
            if len(nz) == 0: continue
            i = r + nz[0]
            if i != r: M[[r, i]] = M[[i, r]]
            M[r] = MUL[s.INV[M[r, c]], M[r]]
            col = M[:, c].copy(); col[r] = 0; rows = np.nonzero(col)[0]
            if len(rows): M[rows] ^= MUL[col[rows][:, None], M[r][None, :]]
            piv.append(c); r += 1
        return M[:r], piv
    def rank(s, M): return len(s.rref(M)[1])
    def nullspace(s, M):
        """rows x with M x = 0"""
        R, piv = s.rref(M); nc = M.shape[1]; sp = set(piv)
        free = [c for c in range(nc) if c not in sp]
        B = np.zeros((len(free), nc), dtype=np.int16)
        for t, f in enumerate(free):
            B[t, f] = 1
            for i, pc in enumerate(piv): B[t, pc] = R[i, f]   # char 2: -x = x
        return B
    def matmul(s, A, B):
        A = np.asarray(A); B = np.asarray(B)
        out = np.zeros((A.shape[0], B.shape[1]), dtype=np.int16)
        for t in range(A.shape[1]): out ^= s.MUL[A[:, t][:, None], B[t][None, :]]
        return out

def model(e, j, m0=0, verbose=True):
    q = 2 ** e; k = q // 2; N = q * q + q + 1
    F = GF(2, e); m, poly = PRIM[e]; B = BF(m, poly)
    # embed F_q into B
    irr = IRR[(2, e)]
    def ev(th, cs):
        r = 0; p = 1
        for c in cs:
            if c: r ^= p
            p = int(B.MUL[p, th])
        return r
    th = next(t for t in range(2, B.n) if ev(t, list(irr)) == 0)
    phi = np.zeros(q, dtype=np.int16)
    for a in range(q):
        v = 0; p = 1
        for i in range(e):
            if a >> i & 1: v ^= p
            p = int(B.MUL[p, th])
        phi[a] = v
    for a in range(q):
        for b in range(q): assert phi[F.mul[a, b]] == B.MUL[phi[a], phi[b]]
    phinv = {int(phi[a]): a for a in range(q)}
    # Singer matrix: g = alpha (=2), g^3 = k0 + k1 g + k2 g^2
    g = 2; g2 = int(B.MUL[g, g]); g3 = int(B.MUL[g2, g])
    sol = [(c0, c1, c2) for c0 in range(q) for c1 in range(q) for c2 in range(q)
           if (int(phi[c0]) ^ int(B.MUL[phi[c1], g]) ^ int(B.MUL[phi[c2], g2])) == g3]
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
    Ainv = minv(A); AT = A.T.copy(); AmT = minv(AT)  # A^{-T}
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
    Dp = Dall[1:]
    # global incidence check: p_i in l_j iff (j - i) mod N in Dall
    Ds = set(Dall)
    for i in range(0, N, 37):
        for jj in range(N): assert (dot(P[i], Lr[jj]) == 0) == (((jj - i) % N) in Ds)
    def coords(u, v, y):
        for lam in range(q):
            for mu in range(q):
                if tuple(int(F.add[F.mul[lam, a], F.mul[mu, b]]) for a, b in zip(u, v)) == tuple(y): return lam, mu
        raise ValueError
    def mono(lam, mu): return [int(F.mul[F.pw(lam, k - 1 - r), F.pw(mu, r)]) for r in range(k)]
    fr = lambda x: F.frob(x, j)
    def lperp(l):
        return [x for x in [(1, 0, 0), (0, 1, 0), (0, 0, 1)] + [(a, b, c) for a in range(q) for b in range(q) for c in range(q)]
                if dot(x, l) == 0 and any(x)]
    # point basis at p_i: A^{-T i}(u0,v0); line basis at l_j: A^j(u0',v0')
    u0, v0 = (0, 1, 0), (0, 0, 1)
    cand = lperp(l0); u0p = cand[0]; v0p = next(x for x in cand if norm(x) != norm(u0p))
    Pu = [u0]; Pv = [v0]; Lu = [u0p]; Lv = [v0p]
    for i in range(1, N):
        Pu.append(mv(AmT, Pu[-1])); Pv.append(mv(AmT, Pv[-1])); Lu.append(mv(A, Lu[-1])); Lv.append(mv(A, Lv[-1]))
    # explicit per-edge local-code coefficients
    gE = {}; hE = {}
    for i in range(N):
        for d in Dp:
            jj = (i + d) % N
            lam, mu = coords(Pu[i], Pv[i], Lr[jj]); gE[(i, d)] = mono(lam, mu)
            lam, mu = coords(Lu[jj], Lv[jj], P[i]); hE[(i, d)] = mono(fr(lam), fr(mu))
    G0 = {d: gE[(0, d)] for d in Dp}; H0 = {d: hE[(0, d)] for d in Dp}
    # determine cocycle: row (i,d) = (s G0_d, t H0_d)
    theta = None
    def ratio(x, y):
        nzi = next(r for r in range(k) if y[r]); s = int(F.mul[x[nzi], F.inv[y[nzi]]])
        assert all(int(F.mul[s, y[r]]) == x[r] for r in range(k)); return s
    for (i, d) in gE:
        s = ratio(gE[(i, d)], G0[d]); t = ratio(hE[(i, d)], H0[d])
        w = (i + d) >= N
        rt = int(F.mul[t, F.inv[s]])
        if not w: assert rt == 1, (i, d, s, t)
        else:
            if theta is None: theta = rt
            assert rt == theta
    if theta is None: theta = 1
    # rho^N = theta (in B)
    lt = int(B.LOG[phi[theta]]); assert lt % N == 0; rho = int(B.EXP[lt // N])
    Gm = np.array([[phi[x] for x in G0[d]] for d in Dp], dtype=np.int16)
    Hm = np.array([[B.MUL[B.pw(rho, d), phi[x]] for x in H0[d]] for d in Dp], dtype=np.int16)
    zeta = int(B.EXP[q - 1])  # order N
    assert B.pw(zeta, N) == 1 and all(B.pw(zeta, N // pp) != 1 for pp in [3, 7, 13, 73] if N % pp == 0)
    if verbose: print(f"[model] q={q} j={j} m0={m0} N={N} Dp={Dp} theta={theta} A=({k0},{k1},{k2})", flush=True)
    return dict(F=F, B=B, q=q, k=k, N=N, e=e, j=j, Dp=Dp, P=P, L=Lr, gE=gE, hE=hE, G=Gm, H=Hm, zeta=zeta,
                rho=rho, theta=theta, phi=phi, norm=norm)

def blocks(Mo):
    """M_chi = [G | diag(zeta^{chi d}) H] for chi in Z_N; returns list of (K_chi rows, Omega_chi rows)."""
    B = Mo['B']; N = Mo['N']; Dp = Mo['Dp']; G = Mo['G']; H = Mo['H']; z = Mo['zeta']
    out = []
    for chi in range(N):
        lam = np.array([B.pw(z, (chi * d) % N) for d in Dp], dtype=np.int16)
        M = np.hstack([G, B.MUL[lam[:, None], H]])
        K = B.nullspace(M); Om = B.nullspace(M.T.copy())
        assert K.shape[0] == Om.shape[0]
        out.append((K, Om, lam))
    return out
