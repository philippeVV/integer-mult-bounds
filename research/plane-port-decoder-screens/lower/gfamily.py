"""Indicative check of the 'incumbent G-decoder + lost-edge correction' family for given planes."""
import json, sys
from ports import P, dot
from gamma import edges, NB, n, orth
emap = {int(k): v for k, v in json.load(open('e8map.json')).items()}
R = [emap[q] for q in P]
ip = lambda a, b: sum(x * y for x, y in zip(a, b)) // 4
PR = 1000003
def rank_mod(M):
    M = [[x % PR for x in row] for row in M]; rk = 0; cols = len(M[0])
    for c in range(cols):
        piv = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]; inv = pow(M[rk][c], PR - 2, PR)
        M[rk] = [x * inv % PR for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]; M[i] = [(x - f * y) % PR for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk
def span_rank(vecs):
    """F2 rank of the dot form restricted to span(vecs)."""
    basis = []
    for v in vecs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    G = [[dot(a, b) for b in basis] for a in basis]
    # F2 rank of G
    rows = [int(''.join(map(str, r)), 2) if r else 0 for r in G]; rk = 0; bs = []
    for v in rows:
        for b in bs: v = min(v, v ^ b)
        if v: bs.append(v)
    return len(bs), len(basis)
def run(r):
    E = set(edges(r))
    lost = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in NB[i]:
            a, b = min(i, j), max(i, j)
            if (a, b) not in E: lost[i][j] = ip(R[i], R[j])
    rl = rank_mod(lost)
    # coordinate decoder supports in the frame basis: frame = 8 missing lines
    MISSING = sorted(set(int(k) for k in emap) - set(P))
    frame = [emap[x] for x in MISSING]
    gcost = []
    for f in frame:
        D = [i for i in range(n) if ip(R[i], f) != 0]
        gcost.append((len(D), span_rank([P[i] for i in D] + [r[i] for i in D])[0], span_rank([P[i] for i in D])[0]))
    return len(E), rl, gcost
if __name__ == '__main__':
    for run_ in json.load(open(sys.argv[1])) if sys.argv[1].endswith('runs.json') else [json.load(open(sys.argv[1]))]:
        e, rl, g = run(run_['r'])
        print('Gamma_E edges', e, 'rank_Q(A_lost) >=', rl, 'G-terms (support, plane cost, line cost):', g, flush=True)
