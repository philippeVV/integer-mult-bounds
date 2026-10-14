"""Exact (Fraction / GF(2)) re-implementation of cold44x.py's per-block pairing computation, instrumented.

Same definitions as research/gauge44-leaver8/construction/notri/cold44x.py (lines building M, MI, cur, protected, I, C,
L, A, K, BK), returned as data instead of constructing the word.  Read-only use of the repo.
"""
import gzip, json, copy
from fractions import Fraction as Q
from pathlib import Path


def rd(p):
    p = Path(p)
    return json.loads(gzip.decompress(p.read_bytes()) if p.suffix == '.gz' else p.read_bytes())


def basis(vs):
    d = {}
    for a in vs:
        for k in sorted(d, reverse=True):
            if a >> k & 1: a ^= d[k]
        if a:
            k = a.bit_length() - 1
            for j in d:
                if d[j] >> k & 1: d[j] ^= a
            d[k] = a
    return tuple(d[k] for k in sorted(d, reverse=True))


def inside(a, b): return basis(list(a) + list(b)) == tuple(b)


def regs(e): return e[1:3] if e[0] == 'add' else e[1:2]


def inverse(e):
    e = copy.deepcopy(e)
    if e[0] == 'add': e[3] = -e[3]
    else: e[2], e[3] = e[3], e[2]
    return e


def apply(e, M, mapping=None):
    t = mapping[e[1]] if mapping is not None else e[1]
    if e[0] == 'add':
        s = mapping[e[2]] if mapping is not None else e[2]; q = Q(e[3], e[4])
        for j, a in list(M[s].items()):
            x = M[t].get(j, Q(0)) + q * a
            if x: M[t][j] = x
            else: M[t].pop(j, None)
    else:
        q = Q(e[2], e[3]); M[t] = {j: q * a for j, a in M[t].items()}


def rref(A):
    a = [list(map(Q, row)) for row in A]; nr = len(a); nc = len(a[0]) if nr else 0; row = 0; piv = []
    for c in range(nc):
        k = next((i for i in range(row, nr) if a[i][c]), None)
        if k is None: continue
        a[row], a[k] = a[k], a[row]; q = a[row][c]; a[row] = [x / q for x in a[row]]
        for i in range(nr):
            if i != row and a[i][c]: q = a[i][c]; a[i] = [x - q * y for x, y in zip(a[i], a[row])]
        piv.append(c); row += 1
        if row == nr: break
    return a, piv


def rank(A): return len(rref(A)[1]) if A and A[0] else 0


def null(A, nc=None):
    if A: a, piv = rref(A); nc = len(a[0])
    else: a = []; piv = []; assert nc is not None
    free = [j for j in range(nc) if j not in piv]; out = [[Q(0) for j in free] for i in range(nc)]
    for k, c in enumerate(free): out[c][k] = Q(1)
    for i, c in enumerate(piv):
        for k, j in enumerate(free): out[c][k] = -a[i][j]
    return out


def mm(A, B):
    if not A: return []
    if not B: return [[] for row in A]
    bt = list(zip(*B)); return [[sum((a * b for a, b in zip(row, col)), Q(0)) for col in bt] for row in A]


def tr(A): return list(map(list, zip(*A)))


class Word:
    """Forward non-target map P (M), its inverse (MI), final frames (cur) as cold44x.py computes them."""

    def __init__(self, w, hamming):
        self.w = w; v = self.v = w['report']['v']; n = self.n = 2 * v + w['report']['R']
        N = self.N = list(range(v)) + list(range(2 * v, n)); ix = self.ix = {x: i for i, x in enumerate(N)}
        inherited_tail = self.inherited_tail = [w['cleanup'][i][:-1] + [list(f)] for i, f in w.get('early_cleanup', [])]
        mainN = [e for e in w['A'] + w['B'] if not v <= e[1] < 2 * v]
        M = [{i: Q(1)} for i in range(len(N))]
        for e in mainN + inherited_tail: apply(e, M, ix)
        MI = [{i: Q(1)} for i in range(len(N))]
        for e in reversed(mainN + inherited_tail): apply(inverse(e), MI, ix)
        self.M, self.MI = M, MI
        cur = {x: (w['ports'][x],) if x < v else basis(w.get('start_extra', {}).get(str(x), [])) for x in N}
        last = {}
        for k, e in enumerate(w['A'] + w['B'] + inherited_tail):
            for x in regs(e):
                if x in cur: cur[x] = basis(e[-1]); last[x] = k
        self.cur, self.last = cur, last
        self.enter = set(map(int, w.get('start_extra', {})))
        self.protected = [ix[x] for x in N if x < v or x in self.enter]
        self.mods = rd(Path(hamming) / 'PACKET-PROGRAMS.json.gz'); self.proto = rd(Path(hamming) / 'PROTOTYPE.json')

    def block_frame(self, b):
        ss = self.mods[b]['receipt']['source_lines']
        return basis([self.proto['binary_ports'][s] for s in ss])

    def pairing(self, b, I=None, protected=None):
        """Return dict with I, L, K, BK, rank for block b (I: list of register ids; default = cold44x's I)."""
        M, MI, ix, N = self.M, self.MI, self.ix, self.N
        F = self.block_frame(b)
        if I is None: I = [x for x in N if inside(self.cur[x], F)]
        Ii = [ix[x] for x in I]
        prot = self.protected if protected is None else protected
        C = [[M[i].get(j, Q(0)) for j in prot] for i in Ii]
        L = tr(null(tr(C), len(Ii))) if Ii else []
        A = [[MI[j].get(i, Q(0)) for i in Ii] for j in prot]; A = [row for row in A if any(row)]
        K = null(A, len(Ii)) if Ii else []
        BK = mm(L, K) if L and K and K[0] else []
        r = rank(BK) if BK else 0
        return dict(block=b, F=F, I=I, rankC=rank(C) if C else 0, L=L, K=K, BK=BK, r=r,
                    dimL=len(L), dimK=(len(K[0]) if K and K[0] else 0))
