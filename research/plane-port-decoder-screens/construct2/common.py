import sys, os, json, numpy as np
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(PKG, 'gain'))
from base import P, dot, rchoices, n, verify, term_cost, pair_rank
P = list(P)
D = np.array([[bin(a & b).count('1') & 1 for b in range(512)] for a in range(512)], dtype=np.int8)
Pa = np.array(P)
RC = [rchoices(q) for q in P]  # 128 each
def orthmask(r):
    r = np.asarray(r)
    M = (D[Pa][:, Pa] == 0) & (D[Pa][:, r] == 0) & (D[r][:, Pa] == 0) & (D[r][:, r] == 0)
    np.fill_diagonal(M, False)
    return M  # True = free (plane-orth)
def cand_masks(a, r):
    """for port a: matrix [128 x n] of free(a,b) for each candidate r_a."""
    r = np.asarray(r); C = np.array(RC[a])
    base = (D[P[a]][Pa] == 0) & (D[P[a]][r] == 0)  # q_a.q_b, q_a.r_b
    M = base[None, :] & (D[C][:, Pa] == 0) & (D[C][:, r] == 0)
    M[:, a] = False
    return M
