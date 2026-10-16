"""expQT: P112 admits q_b = p_sigma(b) with sigma a perfect matching of the orthogonality graph (q_b . p_b = 0, q.q = 1)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from e8 import frames, LAB, dotq, NL
F = set(frames()[0]); P = [a for a in range(NL) if a not in F]
G = nx.Graph(); G.add_nodes_from(P)
G.add_edges_from((a, b) for i, a in enumerate(P) for b in P[i+1:] if dotq(LAB[a], LAB[b]) == 0)
M = nx.max_weight_matching(G, maxcardinality=True)
sig = {}
for a, b in M: sig[a] = b; sig[b] = a
assert len(sig) == 112 and all(dotq(LAB[a], LAB[sig[a]]) == 0 and dotq(LAB[sig[a]], LAB[sig[a]]) == 1 for a in P)
print('P112 orthogonal perfect matching found:', len(M), 'pairs; deg', min(dict(G.degree()).values()))
print('PASS_EXPQT_MATCH')
