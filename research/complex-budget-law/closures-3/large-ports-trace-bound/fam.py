"""expNM: algebraic large-port families (field-line labels over F_q, q = 2^u).

A port is a projective point x of PG(n-1, q) with B(x,x) != 0.  Its label is the F2-subspace
E_x = F_q * x  (dim u) of F_q^n = F2^{un}, with F2 form <v,w> = Tr_{F_q/F2} B(v,w).
B is either the dot form sum x_i y_i (symmetric, char 2) or the Hermitian form sum x_i y_i^Q (q = Q^2).
Pairing rank PR(x,y) = u if B(x,y) != 0 else 0 (field-line alphabet, expM9 L1').  Window h = u n.
Also: 'plus1' variant = affine points (x,1) of PG(n,q) i.e. E8-style '+1' coordinate (subset of PG(n,q)).
"""
import itertools


class GF:
    POLY = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011}

    def __init__(self, u):
        self.u, self.q = u, 1 << u
        p = self.POLY[u]
        self.exp = [0] * (2 * self.q); self.log = [0] * self.q
        x = 1
        for i in range(self.q - 1):
            self.exp[i] = x; self.log[x] = i
            x <<= 1
            if x & self.q: x ^= p
        for i in range(self.q - 1, 2 * self.q): self.exp[i] = self.exp[i - (self.q - 1)]
        assert len(set(self.exp[:self.q - 1])) == self.q - 1, 'poly not primitive'

    def mul(self, a, b):
        if a == 0 or b == 0: return 0
        return self.exp[self.log[a] + self.log[b]]

    def pw(self, a, e):
        if a == 0: return 0
        return self.exp[(self.log[a] * e) % (self.q - 1)]

    def inv(self, a): return self.exp[(self.q - 1 - self.log[a]) % (self.q - 1)]

    def tr(self, a):
        t, x = 0, a
        for _ in range(self.u):
            t ^= x; x = self.mul(x, x)
        assert t in (0, 1)
        return t


def points(F, n):
    """normalised projective points of PG(n-1,q): first nonzero coordinate = 1."""
    out = []
    for v in itertools.product(range(F.q), repeat=n):
        nz = [c for c in v if c]
        if nz and nz[0] == 1: out.append(v)
    return out


def form(F, kind):
    if kind == 'dot':
        def B(x, y):
            s = 0
            for a, b in zip(x, y): s ^= F.mul(a, b)
            return s
    elif kind == 'herm':
        Q = 1 << (F.u // 2); assert Q * Q == F.q
        def B(x, y):
            s = 0
            for a, b in zip(x, y): s ^= F.mul(a, F.pw(b, Q))
            return s
    else:
        raise ValueError(kind)
    return B


def build(u, n, kind, plus1=False):
    F = GF(u); B = form(F, kind)
    P = points(F, n)
    if plus1:   # affine chart x_n = 1 (E8-style '+1' coordinate)
        P = [p for p in P if p[-1] != 0]
        P = [tuple(F.mul(c, F.inv(p[-1])) for c in p) for p in P]
    P = [p for p in P if B(p, p) != 0]
    return F, B, P


def f2_label(F, x):
    """F2 basis of F_q * x inside F2^{un}: rows = bitvectors; coordinate block i holds x_i * g."""
    rows = []
    for j in range(F.u):
        g = 1 << j
        v = 0
        for i, c in enumerate(x): v |= F.mul(g, c) << (F.u * i)
        rows.append(v)
    return rows


def f2_form(F, n, kind):
    """Gram matrix of <v,w> = Tr B(v,w) on the standard F2 basis of F_q^n (bit index u*i + j <-> e_i * 2^j)."""
    B = form(F, kind); h = F.u * n
    basis = []
    for i in range(n):
        for j in range(F.u):
            basis.append(tuple((1 << j) if k == i else 0 for k in range(n)))
    return [[F.tr(B(a, b)) for b in basis] for a in basis]
