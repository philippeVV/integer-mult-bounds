#!/usr/bin/env python3
"""Literal check of the first-pair trajectory for every (a,s,t) in the catalogue,
using formulas (1),(3) of CONSERVED-QQ16.md with the actual 8x112 R and random
noncommuting port-diagonal permutation routes D_u (u=0 identity) on a finite
regular-representation model of size m. Exact integer/Fraction arithmetic.
Usage: trajectory_check.py ARCHIVE"""
import random, sys
from fractions import Fraction as Q
sys.path.insert(0, sys.argv[2])
import pair_catalogue as pc

A, T, U = pc.A, pc.T, pc.U
m = 7
rng = random.Random(20261014)
cols = pc.query_columns(pc.load(sys.argv[1]))
Rinv = {}
for S in (A, T):  # R_S R_S^T = 4I -> R_S^{-1} = R_S^T/4
    Rinv[tuple(S)] = S
perm = lambda: rng.sample(range(m), m)
D = {0: [list(range(m))] * 112}
for u in range(1, 9):
    D[u] = [perm() for _ in range(112)]
Z = lambda: [0] * m
def apply(P, v):  # (D v)[P[i]] = v[i]
    w = Z()
    for i, x in enumerate(v):
        w[P[i]] += x
    return w
def applyinv(P, v):
    return [v[P[i]] for i in range(m)]
def add(*vs):
    return [sum(z) for z in zip(*vs)]
def sc(c, v):
    return [c * x for x in v]
def query(u, f, ports):  # R D_u restricted to ports: -> 8 x m
    out = [Z() for _ in range(8)]
    for p in ports:
        g = apply(D[u][p], f[p])
        for k in range(8):
            if cols[p][k]:
                out[k] = add(out[k], sc(cols[p][k], g))
    return out
def inject(u, S, q):  # D_u^{-1} iota_S R_S^{-1} q, R_S^{-1}=R_S^T/4
    f = {}
    for p in S:
        v = Z()
        for k in range(8):
            if cols[p][k]:
                v = add(v, sc(Q(cols[p][k], 4), q[k]))
        f[p] = applyinv(D[u][p], v)
    return f
def kU(u, q):
    return {p: applyinv(D[u][p], [sum(Q(cols[p][k], 4) * q[k][i] for k in range(8)) for i in range(m)]) for p in U}
ALL = list(range(112))
def forward(a, x, y):
    f = query(a, x, ALL)
    xa, xt, r = inject(a, A, query(a, x, T)), inject(a, T, f), kU(a, f)
    X, Y = {}, {}
    for p in U:
        X[p] = add(x[p], sc(-1, r[p])); Y[p] = add(y[p], X[p])
    for p in A + T:
        X[p] = xa[p] if p in A else xt[p]; Y[p] = add(x[p], y[p], sc(-1, X[p]))
    return X, Y, f
def backward(s, t, x, y):
    ua = inject(t, A, query(t, y, ALL))
    u = {p: (ua[p] if p in A else y[p]) for p in ALL}
    a, b = query(s, x, T), query(s, u, ALL)
    r = kU(s, a)
    xt = inject(s, T, [add(b[k], sc(12, a[k])) for k in range(8)])
    sub = inject(s, T, [add(b[k], sc(11, a[k])) for k in range(8)])
    X, Y = {}, {}
    for p in U:
        X[p] = add(x[p], sc(-1, y[p])); Y[p] = add(y[p], r[p])
    for p in A:
        X[p] = add(x[p], y[p], sc(-1, u[p])); Y[p] = u[p]
    for p in T:
        X[p] = xt[p]; Y[p] = add(u[p], sc(-1, sub[p]))
    return X, Y
def adjoint(u, q):  # B_u^* q : field p = r_p^T q pulled back by D_u(p)^{-1}
    return {p: applyinv(D[u][p], [sum(cols[p][k] * q[k][i] for k in range(8)) for i in range(m)]) for p in U}
count = 0
for s in range(1, 9):
    for t in range(s + 1, 9):
        e = [[int(k == 0 and i == 0) for i in range(m)] for k in range(8)]
        w = adjoint(t, e)
        bw = adjoint(s, query(s, w, U))
        h = {p: add(sc(48, w[p]), sc(-1, bw[p])) for p in U}
        assert not any(any(r) for r in query(s, h, U))
        Bt = query(t, h, U)
        assert any(any(r) for r in Bt)
        for a in range(9):
            Ba = query(a, h, U)
            x, y = {p: Z() for p in ALL}, {p: Z() for p in ALL}
            for p in U:
                x[p] = h[p]
            ia, it = inject(a, A, Ba), inject(t, A, Bt)
            # X_T = P_{a,T} ell_{a,A} P_{t,A} B_t h
            elA = query(a, it, A)
            xt = inject(a, T, elA)
            for p in A:
                x[p] = sc(-1, it[p]); y[p] = sc(-1, x[p])
            for p in T:
                x[p] = xt[p]; y[p] = sc(-1, xt[p])
            assert any(any(x[p]) for p in T), "initial X_T nonzero"
            X, Y, f = forward(a, x, y)
            assert not any(any(r) for r in f), "f_a = 0"
            X, Y = backward(s, t, X, Y)
            assert all(not any(X[p]) for p in ALL), "all X zero"
            assert all(not any(Y[p]) for p in A + T), "selected Y zero"
            assert all(Y[p] == h[p] for p in U), "Y_U = h"
            count += 1
print("PASS_FIRST_PAIR_TRAJECTORY_ALL", count, "pairs (a,s,t); m =", m)
