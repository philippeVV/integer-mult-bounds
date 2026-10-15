"""Orbit bookkeeping for k=7: (1) sanity of the group table EL (each row an F0-preserving isometry
of the line system, rows distinct, identity present, closed under products with 8 random elements);
(2) |Stab| of every non-covering rep; (3) sum over all reps (covering stabs from s2/) of |EL|/|Stab|
equals C(120,7).  Also checks that the covering/non-covering split agrees with a direct cover test."""
import numpy as np, glob, math, random
from model import GM, F0, COSTAR, FULL
EL = np.load("group_elems.npy").astype(np.int64); n = len(EL); assert n == 172032
G = np.abs(np.array(GM))
for i in range(0, n, 1024):
    E = EL[i:i + 1024]
    assert (G[E[:, :, None], E[:, None, :]] == G[None]).all()
    assert np.isin(E[:, F0], F0).all()
rows = {r.tobytes() for r in EL.astype(np.int8)}; assert len(rows) == n
assert np.arange(120, dtype=np.int8).tobytes() in rows
rnd = random.Random(1)
for _ in range(8):
    g = EL[rnd.randrange(n)]
    P = g[EL]                      # g o h
    assert all(r.tobytes() in rows for r in P.astype(np.int8))
print("EL: 172032 distinct F0-preserving |Gram|-isometries, identity present, closed under 8 random right-mults")
z = np.load("reps7.npz")
def dec(c): return [(c >> (7 * (6 - i))) & 127 for i in range(7)]
non = [dec(c) for c in z["non_code"].tolist()]
cov = z["cov_A"].tolist()
def covers(A):
    u = 0
    for a in A: u |= COSTAR[a]
    return u == FULL
assert all(not covers(A) for A in non) and all(covers(A) for A in cov)
def stab(A):
    inA = np.zeros(120, dtype=bool); inA[A] = True
    return int(inA[EL[:, A]].all(axis=1).sum())
sn = [stab(A) for A in non]
sc = np.concatenate([np.load(f)["stab"] for f in sorted(glob.glob("s2/d*.npz"))]).tolist()
assert len(sc) == len(cov)
tot_c = sum(n // s for s in sc); tot_n = sum(n // s for s in sn)
print("covering orbits", len(cov), "sets", tot_c, "; non-covering orbits", len(non), "sets", tot_n)
print("sum", tot_c + tot_n, "C(120,7)", math.comb(120, 7), "MATCH" if tot_c + tot_n == math.comb(120, 7) else "MISMATCH")
