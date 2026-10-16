"""Contract circuit (incidence + rho) for the expNI 784-arc Singer schedule of PG(2,16)-M, register level,
with C/J lanes, plus the GRS transfer T_G mod p for random integer gates (read-only use of ../data/schedule784.json, a copy of expNI's)."""
import json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
P = 33554393  # prime < 2^25: products < 2^50 fit int64

def load(path=os.path.join(HERE, '..', 'data', 'schedule784.json')):   # expNI schedule784.json
    d = json.load(open(path))
    order = d['order']; arcs = [tuple(a) for a in d['arcs']]; R = d['R']
    pos = {v: i for i, v in enumerate(order)}
    return order, arcs, R, pos

class Circuit:
    """Sessions: R arcs minus lanes. mode 'same': C=J=L (12 arcs, rho=id on the other 772).
    mode 'pair': C and J disjoint; R minus (C u J) identity, J-arc sources -> C-arc sinks (paired in given order)."""
    def __init__(self, nlanes=12, seed=0, mode='same'):
        order, arcs, R, pos = load()
        self.order, self.arcs, self.pos = order, arcs, pos
        rng = np.random.default_rng(seed)
        Rl = list(R); perm = rng.permutation(len(Rl))
        if mode == 'same':
            C = set(Rl[i] for i in perm[:nlanes]); J = set(C)
        else:
            C = set(Rl[i] for i in perm[:nlanes]); J = set(Rl[i] for i in perm[nlanes:2 * nlanes])
        self.C, self.J = C, J
        self.src = [r for r in Rl if r not in C]          # column index -> arc whose source feeds it
        self.snk_arcs = [r for r in Rl if r not in J]     # sink arcs (rows before rho-reindex)
        # rho: source column k (arc a) -> sink arc rho(a)
        Jl = sorted(J - C); Cl = sorted(C - J)
        rho = {a: a for a in Rl if a not in C and a not in J}
        for a, b in zip(Jl, Cl): rho[a] = b
        self.rho = rho
        sidx = {a: i for i, a in enumerate(self.snk_arcs)}
        # reindexed rows: row k = sink rho(src[k]) so that the contract identity is diagonal
        self.row_of_sinkarc = {rho[a]: k for k, a in enumerate(self.src)}
        self.W = len(self.src); assert self.W == len(self.snk_arcs)
        self.Rset = set(R)
        V = len(order)
        self.ins = {v: [] for v in range(V)}; self.outs = {v: [] for v in range(V)}
        for i, (u, w) in enumerate(arcs):
            self.outs[u].append(i); self.ins[w].append(i)
        self.colidx = {a: k for k, a in enumerate(self.src)}

    def transfer(self, rng, p=P, gates=None, low=0, high=None):
        """A[k_row, k_col] = T_G[sink rho(src k_row), source k_col]  (rho-reindexed), mod p."""
        W = self.W; val = {}
        A = np.zeros((W, W), dtype=np.int64)
        for v in self.order:
            X = np.zeros((len(self.ins[v]), W), dtype=np.int64)
            for r, a in enumerate(self.ins[v]):
                if a in self.Rset:
                    if a not in self.C: X[r, self.colidx[a]] = 1
                else:
                    X[r] = val.pop(a)
            if gates is not None: G = gates[v]
            else: G = rng.integers(low, p if high is None else high, size=(len(self.outs[v]), len(self.ins[v])), dtype=np.int64)
            Y = matmul(G % p, X, p)
            for r, a in enumerate(self.outs[v]):
                if a in self.Rset:
                    if a not in self.J: A[self.row_of_sinkarc[a]] = Y[r]
                else: val[a] = Y[r]
        assert not val
        return A

def matmul(A, B, p=P):
    out = np.zeros((A.shape[0], B.shape[1]), dtype=np.int64)
    for k in range(A.shape[1]):
        out = (out + np.outer(A[:, k], B[k])) % p
    return out

def rank_mod(M, p=P):
    import flint
    return flint.nmod_mat(M.tolist(), p).rank()

def Zmat(A, Min, Mout, p=P):
    """rho-reindexed: Z[j,k] = A[j,k] (Mout[j] - Min[k]); Min, Mout arrays (W,m,m) mod p."""
    W, m = Min.shape[0], Min.shape[1]
    D = (Mout[:, None, :, :] - Min[None, :, :, :]) % p          # (W,W,m,m)
    Z = (A[:, :, None, None] * D) % p
    return Z.transpose(0, 2, 1, 3).reshape(W * m, W * m)
