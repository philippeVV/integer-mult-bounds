import numpy as np
from xk_sparse import orth
def flatten(q, it=4000):
    d = np.zeros(q.shape[0])
    for _ in range(it):
        Yq = orth(np.exp(d)[:, None] * q); a = np.sum(abs(Yq) ** 2, 1); d += 0.5 * (np.log(.5) - np.log(a))
    return orth(np.exp(d)[:, None] * q)
def XY_codes(V, nP, E, Qs):
    X = np.zeros((E, E)); Y = np.zeros((E, E))
    for i, (es, q) in enumerate(zip(V, Qs)):
        P = q @ q.conj().T; O = 4 * abs(P) ** 2; np.fill_diagonal(O, 0); (X if i < nP else Y)[np.ix_(es, es)] = O
    return X, Y
