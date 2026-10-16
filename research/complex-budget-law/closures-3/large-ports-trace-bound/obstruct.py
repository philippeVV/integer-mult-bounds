"""Field-line obstruction (n = 4): in the free graph of nondegenerate field lines over F_q^4 (any q), the common free
neighbourhood of a free pair {x,y} lies in <x,y>^perp, a nondegenerate F_q-plane, whose free graph is a matching.
So a homomorphism nonorth[S] -> polarity graph needs every edge's common neighbourhood graph to be BIPARTITE.
Report edges whose common neighbourhood contains an odd cycle."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import itertools, numpy as np, networkx as nx
from lattices import SYS
for name, f in [('F4', 'maxsub_F4_4.npy'), ('E6', 'maxsub_E6_4.npy'), ('E8', 'heur_E8_4.npy')]:
    L = SYS[name]()[np.load(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), f))]; A = (L @ L.T) != 0; np.fill_diagonal(A, False); V = len(L)
    bad = 0; tot = 0; ex = None
    for x, y in itertools.combinations(range(V), 2):
        if not A[x, y]: continue
        tot += 1
        C = np.nonzero(A[x] & A[y])[0]
        H = nx.Graph(); H.add_nodes_from(C.tolist())
        H.add_edges_from((a, b) for a, b in itertools.combinations(C.tolist(), 2) if A[a, b])
        if not nx.is_bipartite(H):
            bad += 1
            if ex is None:
                cyc = nx.find_cycle(H)  # any cycle
                ex = (x, y, len(C), H.number_of_edges())
    print(f"{name}-{V}: {bad}/{tot} edges have a non-bipartite common neighbourhood (example edge, |common nbhd|, edges: {ex})"
          f" -> {'NO field-line realisation in PG(3,q) for ANY q' if bad else 'no obstruction from this test'}")
print('PASS_EXPQL_OBSTRUCT')
