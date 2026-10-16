"""expRB / C1 first experiment: Singer-equivariant Tanner codes on PG(2,q) (full, n=q+1, q odd) or
PG(2,q)-M (Singer matching, n=q, q even); dim W = sum_{c in Z_N} dim(A ∩ Λ_c B), Λ_c = diag(ζ_N^{c d}).
Char-0 families are computed mod a large prime P split in Q(ζ_L) (dim mod P >= dim over Q(ζ); equality
for almost all P, checked with two primes). Char-2 values: same formulas over F_{2^r}.
usage: c1.py q [full|minus]"""
import sys, math, random, itertools, json
import numpy as np, galois

def isprime(n):
    return n > 1 and all(n % r for r in range(2, int(n**.5) + 1))
def factor_pp(q):
    for p in range(2, q + 1):
        if q % p == 0:
            e = 0; x = q
            while x % p == 0: x //= p; e += 1
            assert x == 1; return p, e
def split_prime(L, lo):
    P = lo - lo % L + 1
    while not isprime(P): P += L
    return P

# ---------- geometry: PG(2,q) Singer model in F_{q^3} ----------
def geometry(q):
    p, e = factor_pp(q); N = q*q + q + 1
    K = galois.GF(p**(3*e)); g = K.primitive_element
    def tr(x): return x + x**q + x**(q*q)
    # D = {d : Tr(g^{-d}) = 0}, representative element u_d = g^{(-d) mod N}
    D = [d for d in range(N) if tr(g**((-d) % N)) == 0]
    assert len(D) == q + 1
    u = {d: g**((-d) % N) for d in D}
    # basis of T0 over F_q
    Fq = [K(0)] + [g**(N*m) for m in range(q - 1)]  # subfield; g^N generates F_q^*
    b1 = u[D[0]]
    b2 = next(u[d] for d in D[1:])  # distinct projective points -> independent
    coord = {}
    for i, s in enumerate(Fq):
        for j, t in enumerate(Fq):
            coord[int(s*b1 + t*b2)] = (i, j)   # indices into Fq list: 0 -> 0, 1+m -> g^{Nm}
    st = {d: coord[int(u[d])] for d in D}
    return dict(p=p, e=e, N=N, K=K, g=g, D=D, st=st, Fq=Fq)

# index encoding of F_q elements: 0 -> zero, 1+m -> (g^N)^m ; log m in Z_{q-1}

# ---------- local codes, returned as k x n matrices over a field F ----------
def rs_hom(F, nodes, k):
    """homogeneous RS: rows s^{k-1-r} t^r at node (s,t)."""
    return F([[ (s**(k-1-r)) * (t**r) for (s, t) in nodes] for r in range(k)])

def lift_nodes(F, G, Dn, omega, u=1):
    """Teichmüller lift of (s,t) in F_q^2 -> (τ(s)^u, τ(t)^u), τ(g^{Nm}) = ω^m, τ(0)=0."""
    def tau(i):
        return F(0) if i == 0 else omega**(((i - 1) * u) % (G['Fq'].__len__() - 1))
    return [(tau(G['st'][d][0]), tau(G['st'][d][1])) for d in Dn]

def char2_nodes(G, Dn, frob=0):
    """nodes in F_q ⊂ F_{q^3} (char p), Frobenius twist x -> x^{p^frob}."""
    K = G['K']; Fq = G['Fq']; pf = G['p']**frob
    return [(Fq[G['st'][d][0]]**pf, Fq[G['st'][d][1]]**pf) for d in Dn]

def freq_code(F, zeta, Dn, freqs):
    return F([[zeta**((f*d) % len_) for d in Dn] for f in freqs]) if False else \
           F([[zeta**(f*d) for d in Dn] for f in freqs])

# ---------- Tanner dim ----------
def tanner(F, A, B, xi, Dn, N, Q1):
    """xi has order N*Q1; class m in Z_Q1 = edge cocycle rho=xi^m (rho^N in mu_Q1); returns per-class S."""
    kA, kB = A.shape[0], B.shape[0]; M_ = N*Q1; Sm = [0]*Q1; mx = [0]*Q1
    for x in range(M_):
        lam = F([xi**((x*d) % M_) for d in Dn])
        r = np.linalg.matrix_rank(np.vstack([A, B * lam]))
        dd = kA + kB - r; Sm[x % Q1] += dd; mx[x % Q1] = max(mx[x % Q1], dd)
    return Sm, mx

def mds_info(F, A, maxsub=3000, rng=random.Random(1)):
    k, n = A.shape
    subs = list(itertools.combinations(range(n), k))
    if len(subs) > maxsub: subs = rng.sample(subs, maxsub)
    bad = sum(1 for T in subs if np.linalg.matrix_rank(A[:, list(T)]) < k)
    return bad, len(subs), int(np.linalg.matrix_rank(A))

def run(q, model):
    G = geometry(q); N = G['N']; p = G['p']; D = G['D']
    Dn = D[1:] if model == 'minus' else D[:]   # Singer matching: drop D[0]
    n = len(Dn); assert n % 2 == 0; k = n // 2; E = N * n
    out = dict(q=q, model=model, N=N, n=n, k=k, E=E, rows=[])
    Q1 = math.gcd(q - 1, N) if q > 2 else 1  # cocycle classes F_q^*/(F_q^*)^N (rest are gauge-equivalent)
    def rec(name, F, A, B, xi, note=''):
        Sm, mx = tanner(F, A, B, xi, Dn, N, Q1)
        S = max(Sm); m = Sm.index(S)
        bA = mds_info(F, A); bB = mds_info(F, B)
        r = dict(fam=name, char=F.characteristic, S=S, S_by_cocycle=Sm, SE=round(S/E,4), maxblock=mx[m],
                 nonMDS_A=f"{bA[0]}/{bA[1]}", rankA=bA[2], nonMDS_B=f"{bB[0]}/{bB[1]}", rankB=bB[2], note=note)
        out['rows'].append(r); print(json.dumps(r), flush=True)
    # ---- char 0 (mod P) ----
    L = math.lcm(N*Q1, q - 1, p, 8)
    for P in [split_prime(L, 10**9), split_prime(L, 2*10**9 - 10**8)][:2]:
        F = galois.GF(P); h = F.primitive_element
        xi = h**((P - 1)//(N*Q1)); zeta = xi**Q1; omega = h**((P - 1)//(q - 1)); zp = h**((P - 1)//p)
        rng = random.Random(7)
        units = [u for u in range(1, q - 1) if math.gcd(u, q - 1) == 1] or [1]
        for u in units:
            A = rs_hom(F, lift_nodes(F, G, Dn, omega, 1), k)
            B = rs_hom(F, lift_nodes(F, G, Dn, omega, u), k)
            rec(f'teich_RS_twist{u}', F, A, B, xi, f'P={P}')
        for sd in range(2):
            A = F.Random((k, n), seed=sd); B = F.Random((k, n), seed=100 + sd)
            rec(f'generic{sd}', F, A, B, xi, f'P={P}')
        # Singer-frequency codes
        fsets = {'consec': list(range(k)), 'rand': rng.sample(range(N), k)}
        for nm, FA in fsets.items():
            for u in [1, -1] + [x for x in [2, p] if math.gcd(x, N) == 1]:
                FB = [(u*f) % N for f in FA]
                A = freq_code(F, zeta, Dn, FA); B = freq_code(F, zeta, Dn, FB)
                rec(f'freq_{nm}_gal{u}', F, A, B, xi, f'P={P}')
        # extended QR codes on P^1(F_p) (q = p prime), affine coord x = s/t
        if G['e'] == 1 and model == 'full':
            QR = sorted({(x*x) % p for x in range(1, p)}); NQ = [x for x in range(1, p) if x not in QR]
            Fq = G['Fq']; Kf = galois.GF(p)
            def aff(d):
                s, t = (Fq[i] for i in G['st'][d])
                return None if t == 0 else int(Kf(int(s)) / Kf(int(t))) if False else None
            # need s/t in F_p: map subfield elements to integers via log
            lg = {0: None}
            gp = int(galois.GF(p).primitive_element)
            def toint(i):  # index -> integer in F_p: 1+m -> gp^m
                return 0 if i == 0 else pow(gp, i - 1, p)
            xs = []
            for d in Dn:
                si, ti = G['st'][d]
                s, t = toint(si), toint(ti)
                xs.append(None if t == 0 else (s * pow(t, -1, p)) % p)
            for gam in [1, 2, -1]:
                for (lA, lB) in [('Q', 'Q'), ('Q', 'N'), ('N', 'N')]:
                    def qr(sel):
                        Ls = [0] + [(-x) % p for x in (QR if sel == 'Q' else NQ)]
                        rows = []
                        for l in Ls:
                            rows.append([ (p*gam if l == 0 else 0) if x is None else zp**((l*x) % p) for x in xs])
                        return F(np.array([[int(F(v)) if not isinstance(v, int) else v % P for v in r] for r in rows]))
                    rec(f'extQR_{lA}{lB}_gam{gam}', F, qr(lA), qr(lB), xi, f'P={P}')
        if P: pass
    # ---- char 2 / char p controls ----
    if p == 2 or True:
        r = 1
        while pow(2, r, math.lcm(N, q - 1 if p == 2 else 1)) != 1: r += 1
        if p == 2:
            K = G['K']; g = G['g']; xi = g**((q - 1)//Q1)
            for f in range(G['e']):
                A = rs_hom(K, char2_nodes(G, Dn, 0), k); B = rs_hom(K, char2_nodes(G, Dn, f), k)
                rec(f'char2_RS_frob{f}', K, A, B, xi, 'control (expQK formula)')
        # char-2 reductions of the Singer-frequency codes (canonical: N odd, unramified)
        r = 1
        while pow(2, r, N*Q1) != 1: r += 1
        F2 = galois.GF(2**r); x2 = F2.primitive_element**((2**r - 1)//(N*Q1)); z2 = x2**Q1
        rng = random.Random(7)
        fsets = {'consec': list(range(k)), 'rand': rng.sample(range(N), k)}
        for nm, FA in fsets.items():
            for u in [1, -1] + [x for x in [2, p] if math.gcd(x, N) == 1]:
                FB = [(u*f) % N for f in FA]
                rec(f'freq_{nm}_gal{u}_mod2', F2, freq_code(F2, z2, Dn, FA), freq_code(F2, z2, Dn, FB), x2)
    return out

if __name__ == '__main__':
    q = int(sys.argv[1]); model = sys.argv[2] if len(sys.argv) > 2 else ('minus' if q % 2 == 0 else 'full')
    out = run(q, model)
    json.dump(out, open(f'{__file__.rsplit("/",1)[0]}/c1_q{q}_{model}.json', 'w'), indent=1)
    nq = q if q % 2 == 0 else q + 1; Eq = (q * q + q + 1) * nq
    s0 = max((r['S'] for r in out['rows'] if r['char'] != 2 and r['nonMDS_A'].startswith('0/') and r['nonMDS_B'].startswith('0/')), default=0)
    print(f'best char-0 S with MDS gates = {s0}, S/E = {s0 / Eq:.4f} (threshold 1/6)')
    print('PASS_EXPRB_C1' if s0 / Eq < 1 / 6 else 'FAIL_EXPRB_C1')
