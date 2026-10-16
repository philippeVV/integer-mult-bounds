"""Exact verification of a complex-orientation unit (q, w odd in F2^5, decoder N = R A^T over Q) and its bit twin.
complex: N_TS == 0 whenever T != S and w_T.q_S = 1, N_TT != 0 (Fractions); x_c = sum_j rank_F2 B[R_j, A_j] / v, h_c = rank_F2 B.
bit twin (expNB §3 role swap): labels = Q matrix N (h_b = rank_Q N), F2 decoder = B; twin needs even Q-norms (L2)."""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)   # unit.py (expNA)
from unit import rank_f2_rows, rank_q

def dot(a, b): return bin(a & b).count('1') & 1

def verify(hit):
    q, w, R, A = hit['q'], hit['w'], hit['R'], hit['A']
    v = len(q); k = len(R[0])
    assert all(bin(x).count('1') & 1 for x in q + w), 'non-odd label/sink'
    B = [[dot(w[t], q[s]) for s in range(v)] for t in range(v)]
    N = [[sum(Fr(R[t][j]) * A[s][j] for j in range(k)) for s in range(v)] for t in range(v)]
    for t in range(v):
        if B[t][t] != 1 or N[t][t] == 0: return dict(ok=False, why='diag at %d' % t)
        for s in range(v):
            if s != t and B[t][s] and N[t][s] != 0: return dict(ok=False, why='lossy %d<-%d' % (t, s))
    costs = []
    for j in range(k):
        Rj = [t for t in range(v) if R[t][j]]; Aj = [s for s in range(v) if A[s][j]]
        rows = [sum(B[t][s] << i for i, s in enumerate(Aj)) for t in Rj]
        costs.append(rank_f2_rows(rows) if Rj and Aj else 0)
    hc = rank_f2_rows([sum(B[t][s] << s for s in range(v)) for t in range(v)])
    hb = rank_q(N)
    return dict(ok=True, v=v, h_c=hc, cost=sum(costs), costs=costs, x_c=Fr(sum(costs), v), h_b=hb,
                diagN=[str(N[t][t]) for t in range(v)])

if __name__ == '__main__':
    import json
    f = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'hit_v6_k3_c1.json')
    r = verify(json.load(open(f)))
    print(r)
    # sanity hit from cpunit.py with x <= 1 allowed: exact unit at x = 1 (no unit with x < 1 is known)
    assert r['ok'] and r['x_c'] == 1 and r['h_c'] == 4 and r['h_b'] == 3
    print('PASS_EXPNF_VERIFY')
