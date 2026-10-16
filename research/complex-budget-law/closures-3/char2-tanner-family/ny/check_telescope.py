"""Control: verify the telescoping identity  sum_e (b_e a_e^T) (x) (M_head - M_tail) = Z_G
on scenario A with random integer gates and random (non-commuting) frames at every gate and terminal, mod p.
Also reports generic rank of T_G and the number of gate parameters vs. reachable (sink, source) pairs."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import numpy as np
from grs import *

rng = np.random.default_rng(7)
m = 2
gates = {v: rng.integers(0, P, size=(len(PORTS[v][1]), len(PORTS[v][0])), dtype=np.int64) for v in range(182)}
T = transfer(rng, gates=gates)
# forward a, backward b
a = {}
for v in bh.PI:
    ins, outs = PORTS[v]
    X = np.zeros((len(ins), NB), dtype=np.int64)
    for r, q in enumerate(ins):
        X[r] = a[q] if q in PHYS else np.eye(NB, dtype=np.int64)[IDX[q]]
        if q not in PHYS: a[('src', q)] = X[r]
    Y = matmul(gates[v], X)
    for r, q in enumerate(outs):
        if q in PHYS: a[q] = Y[r]
        else: a[('snk', q)] = Y[r]
b = {}
for v in reversed(bh.PI):
    ins, outs = PORTS[v]
    Bo = np.zeros((len(outs), NB), dtype=np.int64)
    for r, q in enumerate(outs):
        Bo[r] = b[q] if q in PHYS else np.eye(NB, dtype=np.int64)[IDX[q]]
        if q not in PHYS: b[('snk', q)] = Bo[r]
    Bi = matmul(gates[v].T.copy(), Bo)
    for r, q in enumerate(ins):
        if q in PHYS: b[q] = Bi[r]
        else: b[('src', q)] = Bi[r]
Fg = {v: rng.integers(0, P, size=(m, m), dtype=np.int64) for v in range(182)}
Min = [rng.integers(0, P, size=(m, m), dtype=np.int64) for _ in range(NB)]
Mout = [(Min[k] + np.eye(m, dtype=np.int64)) % P for k in range(NB)]
Mout = [Mout[k] for k in range(NB)]
LHS = np.zeros((NB * m, NB * m), dtype=np.int64)
def add(bv, av, dM):
    global LHS
    outer = np.zeros((NB, NB), dtype=np.int64)
    for j in np.nonzero(bv)[0]:
        outer[j] = mulmod(np.int64(bv[j]), av)
    for j in np.nonzero(bv)[0]:
        for k in np.nonzero(outer[j])[0]:
            LHS[j*m:(j+1)*m, k*m:(k+1)*m] = (LHS[j*m:(j+1)*m, k*m:(k+1)*m] + mulmod(np.int64(outer[j, k]), dM % P)) % P
for q in bh.phys:
    t, h = bh.tail_head(q)
    add(b[q], a[q], Fg[h] - Fg[t])
for q in B:
    k = IDX[q]
    add(b[('src', q)], a[('src', q)], Fg[bh.sess_src[q]] - Min[k])
    add(b[('snk', q)], a[('snk', q)], Mout[k] - Fg[bh.sess_snk[q]])
Z = Zmat(T, Min, Mout)
ok = np.array_equal(LHS % P, Z % P)
print("telescoping identity holds mod p:", ok)
print("rank T_G mod p =", rank_mod(T), "of", NB, "; nonzero entries of T_G:", int((T != 0).sum()),
      "; gate parameters:", NPAR)
print("rank Z (random generic frames, m=2):", rank_mod(Z), "of", NB * m)
print('PASS_EXPNY_TELESCOPE' if ok else 'FAIL_EXPNY_TELESCOPE')
