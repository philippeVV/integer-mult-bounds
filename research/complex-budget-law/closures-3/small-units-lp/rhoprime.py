"""Exact idempotent-rank cap rho'_K of a block (W, Q): the largest rank of an idempotent X = D_r Bt[R, A] D_a whose
rows/cols use sinks in W and labels in Q.  For row support W' and col support Q' (r, a nonzero there), X^2 = X iff
  Bt[w, q'] = sum_{s=(q_s,w_s): q_s in Q', w_s in W'} lam_s Bt[w, q_s] Bt[w_s, q']   for all w in W', q' in Q'
(lam_s = sum of r a over ports of type s), and then rank X = rank Bt[W', Q'].  Exact over Q with FLINT (nmod check first).
We scan supports W' subset W, Q' subset Q with at most `rem` removed rows/cols, in decreasing rank order.
usage: rhoprime.py n orbit [rem]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, time, itertools, core6, flint
import numpy as np
dot = core6.dot
P = 2**31 - 1
def rank_mod(M):
    return flint.nmod_mat(M, P).rank() if M else 0
def system(W, Q):
    S = [(q, w) for w in W for q in Q if dot(w, q)]
    A = [[dot(w, s[0]) * dot(s[1], qq) for s in S] for w in W for qq in Q]
    b = [dot(w, qq) for w in W for qq in Q]
    return S, A, b
def feasible(W, Q, exact=True):
    S, A, b = system(W, Q)
    if not S: return False
    r1 = flint.nmod_mat(A, P).rank(); r2 = flint.nmod_mat([r + [x] for r, x in zip(A, b)], P).rank()
    if not exact: return r1 == r2
    # mod-p consistent is necessary for Q-consistency only up to bad primes; confirm exactly over Q both ways
    q1 = flint.fmpq_mat(A).rank(); q2 = flint.fmpq_mat([r + [x] for r, x in zip(A, b)]).rank()
    return q1 == q2
def rkQ(W, Q): return flint.fmpq_mat([[dot(w, q) for q in Q] for w in W]).rank() if W and Q else 0
def solve(W, Q):
    S, A, b = system(W, Q)
    M = flint.fmpq_mat(A); x = M.solve(flint.fmpq_mat([[v] for v in b])) if False else None
    return S

if __name__ == '__main__':
    n = int(sys.argv[1]); k = int(sys.argv[2]); rem = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    D = core6.load(n); korb = sorted(set(D['corb']))
    i = D['corb'].index(korb[k]); Wm, Qm = D['cols'][i]; W = core6.bits(Wm); Q = core6.bits(Qm)
    full = rkQ(W, Q); t0 = time.time()
    print('orbit %d beta=%d |W|=%d |Q|=%d rank_Q=%d' % (k, D['beta'][i], len(W), len(Q), full), flush=True)
    best = {}; cnt = 0
    for a in range(rem + 1):
        for b in range(rem + 1 - a if False else rem + 1):
            for dw in itertools.combinations(W, a):
                W2 = [w for w in W if w not in dw]
                for dq in itertools.combinations(Q, b):
                    Q2 = [q for q in Q if q not in dq]
                    r = rkQ(W2, Q2)
                    if r <= max(best, default=-1): continue
                    cnt += 1
                    if feasible(W2, Q2, exact=False) and feasible(W2, Q2):
                        best[r] = (W2, Q2); print('  feasible idempotent support: rank %d, removed W %s Q %s' % (r, dw, dq), flush=True)
    print('max idempotent rank found with <= %d removals per side: %s  (%d supports checked, %.0f s)' % (rem, max(best, default=None), cnt, time.time() - t0))
