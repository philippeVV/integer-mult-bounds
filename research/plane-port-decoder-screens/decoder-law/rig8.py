"""Rank-8 masked completions on P120: are they torus orbits of G/2?  (N_ab N_ba = 1/4 on all non-orth pairs
and N = D (G/2) D^-1 for a diagonal D.)  Also: label span of supports of all ray types."""
import numpy as np, sys
from lmg import *
from catalog import OMEGA, SIMPLE
idx, n, REQ = setup(FULL)
G = np.array([[IP[a][b] for b in range(NL)] for a in range(NL)], float)
for sym in (True, False):
    for s in range(3):
        r = lm(FULL, np.ones((n, 8), bool), np.ones((n, 8), bool), seed=s, tries=1, maxfev=600, sym=sym)
        res, V, W = r; N = V @ W.T
        odd = ~REQ
        pr = (N * N.T)[odd]
        # torus test: D = diag; N_ab / (G_ab/2) = d_a/d_b  -> log|ratio| additive on all odd pairs
        ok = np.all(np.abs(np.abs(pr) - 0.25) < 1e-6)
        sgn = np.sign((N / np.where(G == 0, 1, G / 2)))
        print(f'sym={sym} seed {s} res {res:.1e} rank {np.linalg.matrix_rank(N, 1e-7)}  all odd pairs N_ab*N_ba=1/4: {ok}', flush=True)
# label span / rank check of supports for every ray type (fundamental weights) and random dual pairs
def lab_rank(mask): return rank2([LAB[a] for a in bits(mask)])
for i, w in enumerate(OMEGA):
    S_ = supp(w, FULL); print('omega', i, 'support', bin(S_).count('1'), 'label span dim', lab_rank(S_), 'rho', rho(S_, S_))
