"""line systems (one vector per line, integer coords) for expQL."""
import itertools, numpy as np
def canon(vs):
    out, seen = [], set()
    for v in vs:
        v = tuple(v); w = tuple(-x for x in v)
        k = max(v, w)
        if k not in seen: seen.add(k); out.append(k)
    return np.array(out, dtype=np.int64)
def roots_D(n):
    vs = []
    for i, j in itertools.combinations(range(n), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0] * n; v[i] = 2 * si; v[j] = 2 * sj; vs.append(v)
    return vs
def E8():
    vs = roots_D(8) + [list(s) for s in itertools.product((1, -1), repeat=8) if s.count(-1) % 2 == 0]
    return canon(vs)
def E7():
    L = E8(); s = np.ones(8, dtype=np.int64); return L[(L @ s) == 0]          # roots orthogonal to a root -> E7
def E6():
    L = E8(); a = np.ones(8, dtype=np.int64); b = np.array([2,2,0,0,0,0,0,0])
    return L[((L @ a) == 0) & ((L @ b) == 0)]
def F4():
    vs = roots_D(4) + [[2 * (k == i) * s for k in range(4)] for i in range(4) for s in (1, -1)] + [list(s) for s in itertools.product((1, -1), repeat=4)]
    return canon(vs)
def BW16():
    # RM(1,4) weight-8 words; vectors (+-2)^2 0^14 and (+-1)^8 on a word with even # of minus signs (norm 8)
    words = set()
    for c in itertools.product((0, 1), repeat=5):
        w = tuple((c[0] + sum(c[k + 1] * ((x >> k) & 1) for k in range(4))) % 2 for x in range(16))
        if sum(w) == 8: words.add(w)
    vs = []
    for i, j in itertools.combinations(range(16), 2):
        for si in (2, -2):
            for sj in (2, -2):
                v = [0] * 16; v[i] = si; v[j] = sj; vs.append(v)
    for w in words:
        sup = [i for i in range(16) if w[i]]
        for s in itertools.product((1, -1), repeat=8):
            if s.count(-1) % 2 == 0:
                v = [0] * 16
                for i, x in zip(sup, s): v[i] = x
                vs.append(v)
    return canon(vs)
SYS = {'E8': E8, 'E7': E7, 'E6': E6, 'F4': F4, 'BW16': BW16, 'D4': lambda: canon(roots_D(4)), 'D6': lambda: canon(roots_D(6))}
if __name__ == '__main__':
    for k, f in SYS.items():
        L = f(); G = L @ L.T; nrm = np.diag(G)
        print(k, L.shape, 'rank', np.linalg.matrix_rank(L), 'norms', sorted(set(nrm.tolist())), 'nonorth deg', ((G != 0).sum(1) - 1).min(), ((G != 0).sum(1) - 1).max())
