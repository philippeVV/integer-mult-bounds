#!/usr/bin/env python3
"""Finite-model kernel census of the whole retained replica (generous contract of controller.py):
dim J_k(ker O) = rank[O;J_k] - rank O for each bank k, with no records, with PR2174's pair-0 records
(B0 X_0U, B2 X_0U), and with pair-local records on every bank (bank k>=1: B_{2k+1} X_kU, B_{2k+2} X_kU).
Ranks mod two primes (finite-model evidence, not the actual cover).
usage: QQ16_LIB=1 modelcheck.py ARCHIVE m seed"""
import os, sys, time
os.environ['QQ16_LIB'] = '1'
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import controller as K
from fractions import Fraction as Q
m = K.m
t0 = time.time()
basis = [(k, xy, p, i) for k in range(4) for xy in range(2) for p in K.ALL for i in range(m)]
cols = []
recs = {'pair0': [], 'all': []}
for (k, xy, p, i) in basis:
    banks = [K.zero_bank() for _ in range(4)]
    v = K.Z(); v[i] = Q(1)
    banks[k][xy][p] = v
    obs, _, _ = K.run(banks)
    cols.append(obs)
    X0 = banks[0][0]
    r0 = [x for r in K.query(0, X0, K.U) + K.query(2, X0, K.U) for x in r]
    ra = list(r0)
    for b in (1, 2, 3):
        Xb = banks[b][0]
        ra += [x for r in K.query(2 * b + 1, Xb, K.U) + K.query(2 * b + 2, Xb, K.U) for x in r]
    recs['pair0'].append(r0); recs['all'].append(ra)
print('columns built', len(cols), 'obs length', len(cols[0]), f'{time.time()-t0:.0f}s', flush=True)
nb = len(basis)
def Jrows(k):
    rows = []
    for j, (kk, xy, p, i) in enumerate(basis):
        if kk == k and p in K.A + K.T:
            r = [0] * nb; r[j] = 1; rows.append(r)
    return rows
def tomat(colvecs, P):
    out = np.zeros((len(colvecs[0]), len(colvecs)), dtype=np.int64)
    for j, c in enumerate(colvecs):
        for r, x in enumerate(c):
            if x:
                out[r, j] = (x.numerator % P) * pow(x.denominator % P, P - 2, P) % P
    return out
def rank(M, P):
    M = M.copy() % P
    r = 0
    rows, ncol = M.shape
    for c in range(ncol):
        if r == rows: break
        nzr = np.nonzero(M[r:, c])[0]
        if len(nzr) == 0: continue
        piv = r + nzr[0]
        if piv != r: M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), P - 2, P)
        M[r] = (M[r] * inv) % P
        below = np.nonzero(M[r + 1:, c])[0] + r + 1
        if len(below):
            f = M[below, c].copy()
            M[below] = (M[below] - (f[:, None] * M[r][None, :]) % P) % P
        r += 1
    return r
def dedup(M):
    M = M[np.any(M != 0, axis=1)]
    return np.unique(M, axis=0)
for P in (2147483647, 2147483629):
    O = dedup(tomat(cols, P))
    for name, extra in (('none', None), ('pair0 records (16/owner)', recs['pair0']), ('records on all banks (64/owner)', recs['all'])):
        Om = O if extra is None else dedup(np.vstack([O, tomat(extra, P)]))
        rO = rank(Om, P)
        lost = []
        for k in range(4):
            JM = np.array(Jrows(k), dtype=np.int64)
            lost.append(rank(np.vstack([Om, JM]), P) - rO)
        print(f'p={P} m={m} records={name:32s} rank O={rO}  dim J_k(ker O) per bank={lost}  {time.time()-t0:.0f}s', flush=True)
