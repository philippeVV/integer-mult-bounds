"""expQW flat.py: inflation facts behind Lemma S / Theorem W3 on E8 sub-port-sets (P120, P112, P104).

For each port a (F2 label p_a in F2^9, standard dot form) and a Y-label system q_a = p_sigma(a) (sigma an orthogonal
perfect matching), report:
  dense : span{p_b : b free neighbour of a in the set} == p_a^perp  (G/2's transfer row is every free pair, so a
          lossless correction set equal to G/2's transfer row saturates Y_a to p_a^perp: Lemma I)
  flat  : max over hyperplanes K of p_a^perp with q_a in K of #{free neighbours b with p_b in K}
          = the most transfers row a can keep while Y_a stays unsaturated (U_a = q_a + span p_T(a) inside K).
Run: $S/expO/run flat.py   (needs networkx: uv run --with networkx)
"""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from e8 import frames, LAB, dotq, NL, IP, FULL, bits


def rank2(vs):
    basis = []
    for r in vs:
        for b in basis:
            r = min(r, r ^ b)
        if r:
            basis.append(r)
    return len(basis)


def perp(v):
    return [w for w in range(1, 512) if dotq(v, w) == 0]


fr = frames()
f0, f1 = set(fr[0]), None
for f in fr:
    if not (set(f) & f0):
        f1 = set(f); break
SETS = {'P120': list(range(NL)), 'P112': [a for a in range(NL) if a not in f0],
        'P104': [a for a in range(NL) if a not in f0 and a not in f1]}

ok = True
for name, P in SETS.items():
    G = nx.Graph(); G.add_nodes_from(P)
    G.add_edges_from((a, b) for i, a in enumerate(P) for b in P[i + 1:] if dotq(LAB[a], LAB[b]) == 0)
    M = nx.max_weight_matching(G, maxcardinality=True)
    sig = {}
    for a, b in M:
        sig[a] = b; sig[b] = a
    perfect = len(sig) == len(P)
    dense, flat = 0, []
    for a in P:
        pa = LAB[a]
        nb = [b for b in P if b != a and IP[a][b] in (1, -1)]       # free pairs = G/2 transfer row
        assert all(dotq(pa, LAB[b]) == 0 for b in nb)
        sp = rank2([LAB[b] for b in nb])
        dense += (sp == 8)
        if a in sig:
            qa = LAB[sig[a]]
            # hyperplanes K of pa^perp: K = pa^perp cap w^perp for w not in <pa>; distinct K <-> w mod <pa>
            best = 0
            seen = set()
            for w in range(1, 512):
                if w == pa:
                    continue
                key = min(w, w ^ pa)
                if key in seen:
                    continue
                seen.add(key)
                if dotq(qa, w):            # q_a must lie in K
                    continue
                c = sum(1 for b in nb if dotq(LAB[b], w) == 0)
                best = max(best, c)
            flat.append((best, len(nb)))
    print(f'{name}: V={len(P)}  orthogonal perfect matching: {perfect} ({len(M)} pairs)')
    print(f'   dense rows (free-neighbour span = p_a^perp): {dense}/{len(P)}')
    if flat:
        mx = max(f for f, _ in flat); deg = sorted(set(d for _, d in flat))
        print(f'   free degree {deg}; max transfers keeping Y_a unsaturated: {mx} '
              f'(fraction {mx}/{deg[-1]}), distribution {sorted(set(f for f,_ in flat))}')
    ok &= dense == len(P)
print('PASS_EXPQW_FLAT' if ok else 'FAIL')
