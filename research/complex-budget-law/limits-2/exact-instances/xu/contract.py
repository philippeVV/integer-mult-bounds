"""(1) The q = 9 exact schedule as a contract circuit (roles, rho, physical edges), checked against upstream
scripts/finite_bit_contract.py where the binary contract allows it.  Writes contract9.json.
usage: contract.py [path/to/integer-mult-bounds]"""
import os, sys, json
import numpy as np
from collections import deque
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import circ9 as c
REPO = sys.argv[1] if len(sys.argv) > 1 else os.environ.get('IMB_REPO', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), *(['..'] * 5))))
sys.path.insert(0, os.path.join(REPO, 'scripts'))
import finite_bit_contract as fbc

W, V = c.W, c.V
nxt, paths = c.roles()
# role w: source = head of back arc R[w]; its trail; ends at the sink half of back arc R[j]
end_sess = [c.RIDX[end] for _, _, end in paths]
assert sorted(end_sess) == list(range(W))
arc_role = {}
for w, (e0, path, end) in enumerate(paths):
    arc_role[('in', e0)] = w
    for e in path: arc_role[e] = w
    arc_role[('out', end)] = w
def in_role(e): return arc_role[('in', e)] if e in c.RSET else arc_role[e]
def out_role(e): return arc_role[('out', e)] if e in c.RSET else arc_role[e]

gates = []; gate_f9 = []
for v in c.ORDER:
    rin = [in_role(e) for e in c.INS[v]]; rout = [out_role(e) for e in c.OUTS[v]]
    assert rin == rout, 'transition system pairs INS[v][i] with OUTS[v][i]'
    gates.append(dict(roles=rin, xors=[]))
    gate_f9.append(dict(vertex=v, roles=rin, F9_matrix=c.GATE9[v].tolist()))   # out-roles x in-roles, F9 index coding

# ---- upstream checks that do not need XOR gates
fbc.validate_circuit(W, gates)
edges = fbc.physical_graph(W, gates)
n_term_in = sum(1 for u, v, r in edges if u < W); n_term_out = sum(1 for u, v, r in edges if v >= W + len(gates))
print(f"upstream validate_circuit: OK (W = {W}, {len(gates)} gates, every gate on 5 distinct roles)")
print(f"upstream physical_graph: {len(edges)} physical edges = {n_term_in} source-terminal + "
      f"{len(edges) - n_term_in - n_term_out} gate-gate + {n_term_out} sink-terminal "
      f"(expected 166 + 744 + 166 = 1076: {len(edges) == 1076 and n_term_in == n_term_out == 166})")
# gate-gate physical edges are exactly the 744 forward arcs of Gamma - R
gpos = {v: W + i for i, v in enumerate(c.ORDER)}
gg = sorted((u, v) for u, v, r in edges if u >= W and v < W + len(gates))
fw = sorted((gpos[a], gpos[b]) for e, (a, b) in enumerate(c.ARCS) if e not in c.RSET)
print(f"  gate-gate edges == forward arcs of Gamma - R (as a multiset): {gg == fw}")
# first/last gate of each role = head / tail of the role's back arcs
first = {}; last = {}
for u, v, r in edges:
    if u < W: first[r] = v
    if v >= W + len(gates): last[r] = u
ok = all(first[w] == gpos[c.ARCS[paths[w][0]][1]] and last[w] == gpos[c.ARCS[paths[w][2]][0]] for w in range(W))
print(f"  role terminals sit at the back-arc head (source) / tail (sink): {ok}")

# ---- exact F9 scalar map at role level (twin of fbc.scalar_permutation, over F9 instead of XOR/F2)
val = np.full((W, W), c.Z9, dtype=np.int64)
for w in range(W): val[w, w] = c.O9
for g, gf in zip(gates, gate_f9):
    r = g['roles']; val[r] = c.mat_mul9(np.array(gf['F9_matrix']), val[r])
rho = [None] * W; perm_ok = True
for out in range(W):
    nz = np.nonzero(val[out] != c.Z9)[0]
    if len(nz) != 1 or val[out, nz[0]] != c.O9: perm_ok = False; break
    rho[int(nz[0])] = out
print(f"F9 scalar map at role level is a unit-coefficient permutation: {perm_ok}")
expect = [None] * W
for w in range(W): expect[w] = end_sess.index(w)          # session w's value leaves at the sink of role ending at R[w]
print(f"  rho(input role) = output role equals 'role whose trail ends at the same back arc': {rho == expect}")
cyc, seen = [], set()
for w in range(W):
    if w in seen: continue
    L = 0; x = w
    while x not in seen: seen.add(x); x = rho[x]; L += 1
    cyc.append(L)
print(f"  rho fixed points {sum(1 for x in cyc if x == 1)}, cycle type {sorted(cyc, reverse=True)[:12]}...")

# ---- what the binary contract does with this gate list
xor_rho = fbc.scalar_permutation(W, gates)
print(f"upstream scalar_permutation (XOR semantics, no XOR updates expressible): rho = id: {list(xor_rho) == list(range(W))}; "
      f"equals the F9 rho: {list(xor_rho) == rho}")
cand = fbc.zero_frame_candidate(W, gates, m=2); cand['rho'] = rho
try:
    fbc.inspect_candidate(cand, require_deficit=False); print("inspect_candidate accepted (unexpected)")
except ValueError as err:
    print(f"upstream inspect_candidate with the true F9 rho: REJECTED -> {err!r}")

# ---- routing: impossible by girth 6 (any head->tail path in Gamma-R has >= 5 arcs) -- check shortest paths
adj = [[] for _ in range(V)]
for e, (a, b) in enumerate(c.ARCS):
    if e not in c.RSET: adj[a].append(b)
def dist(s, t):
    d = {s: 0}; dq = deque([s])
    while dq:
        x = dq.popleft()
        if x == t: return d[x]
        for y in adj[x]:
            if y not in d: d[y] = d[x] + 1; dq.append(y)
    return None
ds = [dist(c.SRC[w], c.SNK[w]) for w in range(W)]
print(f"routing lower bound: min/max shortest head->tail distance in Gamma-R = {min(ds)}/{max(ds)} arcs, "
      f"sum = {sum(ds)} > 744 forward arcs: {sum(ds) > 744}  => NOT routable (edge-disjoint), girth bound 5*166 = 830")

json.dump(dict(W=W, rho=rho, gates=gates, gate_order=c.ORDER, f9_gates=gate_f9,
               f9_index_note='F9 elements are indices into sorted(Fq) of expXQ/code.py (F_{3^6} model); index 0 = 0, '
                             f'index {c.O9} = 1', role_trails=[[int(e0), [int(x) for x in p], int(end)] for e0, p, end in paths]),
          open(os.path.join(HERE, 'contract9.json'), 'w'))
print("wrote contract9.json")
