"""Snap a numerical in/out code to its support pattern and test (a) the snapped code numerically, (b) the GENERIC dimension of
the support pattern exactly mod two primes (random entries on the support; upper bound on generic char-0 nullity).
dim W = nullity(T - I) with T the back-arc transfer of the DAG Gamma - R.
usage: snap.py NAME solfile tol [tol ...]"""
import sys, json, os, numpy as np
sys.argv, args = sys.argv[:3], sys.argv
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gstruct.py')).read().split("if __name__ == '__main__':")[0])
PR = (1000000007, 998244353)

def rank_mod(M, p):
    M = M.copy() % p; n, c = M.shape; r = 0
    for col in range(c):
        if r == n: break
        nz = np.nonzero(M[r:, col])[0]
        if len(nz) == 0: continue
        piv = r + nz[0]; M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, col]), p - 2, p)) % p
        rows = np.nonzero(M[:, col])[0]; rows = rows[rows != r]
        if len(rows):
            f = M[rows, col][:, None].astype(object); M[rows] = ((M[rows].astype(object) - f * M[r].astype(object)) % p).astype(np.int64)
        r += 1
    return r

def transfer_mod(supp, p, rng):
    Rset = set(R); idx = {e: i for i, e in enumerate(R)}; nR = len(R); val = {}; T = np.zeros((nR, nR), dtype=object); bad = 0
    for v in OR['order']:
        ins, outs, S = supp[v]
        while True:
            Gv = (rng.integers(1, p, size=(k, k)) * S).astype(object)
            if rank_mod(Gv.astype(np.int64), p) == k: break
            bad += 1
            if bad > 200: return None
        X = np.zeros((k, nR), dtype=object)
        for a, e in enumerate(ins):
            if e in Rset: X[a, idx[e]] = 1
            else: X[a] = val[e]
        Y = Gv.dot(X) % p
        for b, e in enumerate(outs):
            if e in Rset: T[idx[e]] = Y[b]
            else: val[e] = Y[b]
    return T

gs0 = gates(Q); out = []
for tol in map(float, args[3:]):
    supp = []; gsn = []
    for ins, outs, Gv in gs0:
        S = np.abs(Gv) > tol * np.abs(Gv).max(); supp.append((ins, outs, S.astype(np.int64))); gsn.append((ins, outs, Gv * S))
    inv_ok = all(abs(np.linalg.det(g)) > 1e-9 for _, _, g in gsn)
    T = transfer(gsn); sv = np.sort(np.linalg.svd(T - np.eye(len(R)), compute_uv=False))
    nnum = int((sv < 1e-9).sum())
    gen = []
    for p in PR:
        Tp = transfer_mod(supp, p, np.random.default_rng(p % 1000))
        gen.append(None if Tp is None else len(R) - rank_mod(((Tp - np.eye(len(R), dtype=np.int64)) % p).astype(np.int64), p))
    rec = dict(sol=sol, tol=tol, snapped_invertible=inv_ok, snapped_nullity_1e9=nnum, sv_edge=[float(x) for x in sv[max(0, nnum - 1):nnum + 2]],
               generic_pattern_nullity_mod_p=gen, nnz=int(sum(s[2].sum() for s in supp)), R=len(R), E_over_6=E / 6)
    print(json.dumps(rec), flush=True)
