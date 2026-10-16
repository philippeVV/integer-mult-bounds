"""expXA common: GF(2) rank, label pairing, masked-decoder validity and cost (expNZ model).

Model (expNZ / expQP): ports a with integral Gram W (W_aa even), labels B = (W + J) mod 2.
A masked decoder is N = sum_k p_k q_k^T over Q with N_aa = 1 and N_ab = 0 whenever B_ab = 1 (a != b).
Cost C = sum_k rank_F2 B[R_k, S_k], R_k = supp p_k, S_k = supp q_k.  gamma = C / V.  h = rank_F2 B.
"""
from fractions import Fraction


def rank_f2_rows(rows):
    """rank over F2 of a list of python-int bit rows."""
    piv = {}
    r = 0
    for v in rows:
        while v:
            lb = v.bit_length() - 1
            if lb in piv:
                v ^= piv[lb]
            else:
                piv[lb] = v
                r += 1
                break
    return r


def block_rank(Bbits, R, S):
    """rank_F2 of B[R,S]; Bbits[a] = int with bit b set iff B_ab = 1."""
    mask = 0
    for b in S:
        mask |= 1 << b
    return rank_f2_rows([Bbits[a] & mask for a in R])


def label_bits(W):
    V = len(W)
    out = []
    for a in range(V):
        x = 0
        Wa = W[a]
        for b in range(V):
            if (Wa[b] + 1) & 1:
                x |= 1 << b
        out.append(x)
    return out


def check_decoder(W, Bbits, terms):
    """terms: list of (coef Fraction, vec p dict a->int, vec q dict a->int). Exact check N_aa=1, N_ab=0 on required."""
    V = len(W)
    N = [dict() for _ in range(V)]
    for c, p, q in terms:
        for a, pa in p.items():
            row = N[a]
            for b, qb in q.items():
                row[b] = row.get(b, 0) + c * pa * qb
    for a in range(V):
        if N[a].get(a, 0) != 1:
            return False, ('diag', a, N[a].get(a, 0))
        for b, v in N[a].items():
            if b != a and v != 0 and (Bbits[a] >> b) & 1:
                return False, ('req', a, b, v)
    return True, N


def cost(Bbits, terms):
    tot = 0
    det = []
    for c, p, q in terms:
        r = block_rank(Bbits, sorted(p), sorted(q))
        det.append(r)
        tot += r
    return tot, det
