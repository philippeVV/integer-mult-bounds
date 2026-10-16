"""Singer difference set for PG(2,q), q prime; normalised so that q*D = D (multiplier fixed)."""
import itertools
def singer(q):
    N = q*q+q+1
    def mul(a, b, f):  # polys as lists len 3 mod monic cubic f=(c0,c1,c2): x^3 = -(c0+c1x+c2x^2)
        r = [0]*5
        for i in range(3):
            for j in range(3): r[i+j] = (r[i+j] + a[i]*b[j]) % q
        for k in (4, 3):
            c = r[k]; r[k] = 0
            for t in range(3): r[k-3+t] = (r[k-3+t] - c*f[t]) % q
        return r[:3]
    for f in itertools.product(range(q), repeat=3):
        if f[0] == 0: continue
        # irreducible: no roots
        if any((x**3 + f[2]*x*x + f[1]*x + f[0]) % q == 0 for x in range(q)): continue
        # order of x must be q^3-1
        p = [1, 0, 0]; x = [0, 1, 0]; pw = []
        for i in range(q**3-1):
            pw.append(p); p = mul(p, x, f)
            if p == [1, 0, 0] and i < q**3-2: break
        else:
            if p != [1, 0, 0]: continue
            D = sorted({i % N for i in range(q**3-1) if pw[i][2] == 0})
            assert len(D) == q+1
            diffs = {(a-b) % N for a in D for b in D if a != b}; assert len(diffs) == N-1
            # normalise: q*D = D + s  -> shift t with s = (q-1) t
            s = next(s for s in range(N) if sorted((q*d) % N for d in D) == sorted((d+s) % N for d in D))
            t = s * pow(q-1, -1, N) % N
            D = sorted((d - t) % N for d in D)
            assert sorted((q*d) % N for d in D) == D
            return N, D
if __name__ == "__main__":
    N, D = singer(17); print(N, D)
    orb = []; seen = set()
    for d in D:
        if d in seen: continue
        o = [d, d*17 % N, d*289 % N]; seen |= set(o); orb.append(o)
    print(orb)

def singer_coords(q):
    """N, D (normalised, qD=D), and for each d in D its projective coordinate (s,t) in F_q^2 on the line
    {a + b x} of F_{q^3}/F_q^* (the same labelling used for point and line stars, as in expRB c1)."""
    import itertools
    N = q*q+q+1
    def mul(a, b, f):
        r = [0]*5
        for i in range(3):
            for j in range(3): r[i+j] = (r[i+j] + a[i]*b[j]) % q
        for kk in (4, 3):
            c = r[kk]; r[kk] = 0
            for t in range(3): r[kk-3+t] = (r[kk-3+t] - c*f[t]) % q
        return r[:3]
    for f in itertools.product(range(q), repeat=3):
        if f[0] == 0 or any((x**3 + f[2]*x*x + f[1]*x + f[0]) % q == 0 for x in range(q)): continue
        p = [1, 0, 0]; pw = []; ok = True
        for i in range(q**3-1):
            pw.append(p); p = mul(p, [0, 1, 0], f)
            if p == [1, 0, 0] and i < q**3-2: ok = False; break
        if not ok: continue
        D0 = sorted({i % N for i in range(q**3-1) if pw[i][2] == 0})
        s = next(s for s in range(N) if sorted((q*d) % N for d in D0) == sorted((d+s) % N for d in D0))
        t = s * pow(q-1, -1, N) % N
        D = sorted((d - t) % N for d in D0)
        st = {}
        for d in D:
            a, b, _ = pw[(d + t) % N]
            # normalise projective (a:b)
            if b: st[d] = ((a * pow(b, -1, q)) % q, 1)
            else: st[d] = (1, 0)
        return N, D, st
