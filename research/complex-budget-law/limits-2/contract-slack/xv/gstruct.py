"""Structure of a found in/out code: gate maps G_v = Q_out Q_in^{-1} (k x k), support at tol, block decomposition,
numerical Tanner dimension, transfer T on the back arcs R and nullity(T - I).
usage: struct.py NAME solfile [tol]"""
import sys, json, os, numpy as np, networkx as nx
HERE = os.path.dirname(os.path.abspath(__file__))
name, sol = sys.argv[1], sys.argv[2]; tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6
G = json.load(open(os.path.join(HERE, '..', 'expXP', 'graphs.json')))[name]
EV = np.array(G['Pts'] + G['Lns']); V, d = EV.shape; k = d // 2; E = G['E']
OR = json.load(open(os.path.join(HERE, f'orient_{name}.json'))); TAIL, HEAD = np.array(OR['tail']), np.array(OR['head']); R = OR['R']
Q = np.load(sol)

def gates(Q):
    out = []
    for v in range(V):
        li = [i for i, e in enumerate(EV[v]) if HEAD[e] == v]; lo = [i for i, e in enumerate(EV[v]) if TAIL[e] == v]
        Gv = Q[v][lo] @ np.linalg.inv(Q[v][li])
        out.append((EV[v][li], EV[v][lo], Gv))
    return out

def analyse(Q, tol=tol, verbose=True):
    gs = gates(Q); blocks = {}; nz = 0
    for ins, outs, Gv in gs:
        sc = np.abs(Gv).max(); S = np.abs(Gv) > tol * sc; nz += S.sum()
        B = nx.Graph(); B.add_nodes_from([('i', a) for a in range(k)] + [('o', b) for b in range(k)])
        for b in range(k):
            for a in range(k):
                if S[b, a]: B.add_edge(('i', a), ('o', b))
        for comp in nx.connected_components(B):
            ni = sum(1 for x in comp if x[0] == 'i'); no = len(comp) - ni
            blocks[(ni, no)] = blocks.get((ni, no), 0) + 1
    return gs, blocks, nz

def tanner_dim(Q, thr=1e-7):
    K = 2 * np.eye(E, dtype=complex)
    for v in range(V):
        K[np.ix_(EV[v], EV[v])] -= Q[v] @ Q[v].conj().T
    w, U = np.linalg.eigh(K); return w, U

def transfer(gs):
    """T: C^R -> C^R through the DAG Gamma - R (source = head of back arc, sink = tail)."""
    Rset = set(R); idx = {e: i for i, e in enumerate(R)}; nR = len(R)
    val = {}  # edge -> row vector over sources
    order = OR['order']; gv = {v: gs[v] for v in range(V)}
    T = np.zeros((nR, nR), dtype=complex)
    for v in order:
        ins, outs, Gv = gv[v]
        X = np.zeros((k, nR), dtype=complex)
        for a, e in enumerate(ins):
            if e in Rset: X[a, idx[e]] = 1
            else: X[a] = val[e]
        Y = Gv @ X
        for b, e in enumerate(outs):
            if e in Rset: T[idx[e]] = Y[b]
            else: val[e] = Y[b]
    return T

if __name__ == '__main__':
    gs, blocks, nz = analyse(Q)
    w, U = tanner_dim(Q)
    T = transfer(gs); sv = np.linalg.svd(T - np.eye(len(R)), compute_uv=False)
    print(json.dumps(dict(sol=sol, tol=tol, blocks={f'{a}x{b}': c for (a, b), c in sorted(blocks.items())}, nnz_gate_entries=int(nz), of=V * k * k,
                          smallest_K=[float(x) for x in w[:3]], K_around=[float(x) for x in w[max(0, int((w < 1e-5).sum()) - 2):int((w < 1e-5).sum()) + 3]],
                          nK_1e5=int((w < 1e-5).sum()), R=len(R), nullity_T_minus_I_1e5=int((sv < 1e-5 * max(1, sv.max())).sum()),
                          sv_TmI_small=[float(x) for x in np.sort(sv)[:3]], sv_TmI_edge=[float(x) for x in np.sort(sv)[int((sv < 1e-4).sum()) - 1:int((sv < 1e-4).sum()) + 2]])))
