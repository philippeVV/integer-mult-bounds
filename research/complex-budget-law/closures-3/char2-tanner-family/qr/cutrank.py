"""expL6 F1 (cut-aligned terminal frames): GRS room per coordinate <= (|A| - rk T[rhoA,A]) + (|Ab| - rk T[rhoAb,Ab]).
Measure the generic rank drop for A = sessions on one side of expQD's near-minimum bisections (several session->side rules),
on schedule784 (+12 C/J lanes, both lane modes) and on expQD's kah1544 / pg1534 order layouts (all back arcs as sessions)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, json, pickle, numpy as np
from circ import Circuit, rank_mod, P, matmul
import os
S_ = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', '')   # expQD layouts, pickles converted to JSON

class GenCircuit(Circuit):
    def __init__(self, order, arcs, R):
        self.order = order; self.arcs = arcs; self.pos = {v: i for i, v in enumerate(order)}
        self.C = set(); self.J = set(); Rl = sorted(R); self.Rset = set(Rl)
        self.src = Rl; self.snk_arcs = Rl; self.rho = {a: a for a in Rl}
        self.row_of_sinkarc = {a: k for k, a in enumerate(Rl)}; self.colidx = dict(self.row_of_sinkarc)
        self.W = len(Rl); V = len(order)
        self.ins = {v: [] for v in range(V)}; self.outs = {v: [] for v in range(V)}
        for i, (u, w) in enumerate(arcs): self.outs[u].append(i); self.ins[w].append(i)

def sides(c, Svec):
    """session k (column) -> side by rule; source port at head(src arc), sink port at tail(rho(src arc))."""
    out = {}
    head = np.array([Svec[c.arcs[a][1]] for a in c.src]); tail = np.array([Svec[c.arcs[c.rho[a]][0]] for a in c.src])
    out['head'] = head == 1; out['tail'] = tail == 1; out['both'] = (head == 1) & (tail == 1)
    return out

def drop(Amat, mask):
    a = np.nonzero(mask)[0]; b = np.nonzero(~mask)[0]
    ra = rank_mod(Amat[np.ix_(a, a)]) if len(a) else 0; rb = rank_mod(Amat[np.ix_(b, b)]) if len(b) else 0
    return len(a), ra, len(b), rb, (len(a) - ra) + (len(b) - rb)

res = []
def run(label, c, cuts, nG=3):
    for g in range(nG):
        A = c.transfer(np.random.default_rng(9000 + g))
        for cname, Svec in cuts.items():
            for rule, mask in sides(c, Svec).items():
                na, ra, nb, rb, d = drop(A, mask)
                rec = dict(layout=label, W=c.W, g=g, cut=cname, rule=rule, A=na, rkA=ra, Abar=nb, rkAbar=rb, drop=d, thresh=round(0.02 * c.W, 1))
                print(json.dumps(rec), flush=True); res.append(rec)

sing = json.load(open(S_ + 'singsa_59_1.json'))
for mode in ['same', 'pair']:
    c = Circuit(mode=mode, seed=0)
    # cyclic-layout halves (Singer orbit union / prefix cut of the 784 order) and expQD's 1542 bisection
    pref = np.zeros(546, dtype=int); pref[c.order[:273]] = 1
    pts = np.zeros(546, dtype=int); pts[:273] = 1   # points vs lines (cut = all 4368 edges: control)
    cuts = {'singsa_59_1(1542)': np.array(sing['S']), 'order_prefix_half': pref, 'points_vs_lines': pts}
    run('schedule784_' + mode, c, cuts)
for f in ['kah1544_ordsa1.json', 'pg1534_ordsa1.json']:
    d = json.load(open(S_ + f))
    if 'order' not in d: print('skip', f, list(d)); continue
    order = list(d['order']); pos = {v: i for i, v in enumerate(order)}; back = set(d['back']); arcs = []
    for ei, (x, y) in enumerate(d['edges']):
        u, w = (x, y) if pos[x] < pos[y] else (y, x)
        arcs.append((w, u) if ei in back else (u, w))
    c = GenCircuit(order, arcs, back)
    cuts = {f.split('_')[0]: np.array(d['S'])}
    run(f.replace('.json', ''), c, cuts)
json.dump(res, open(S_ + '../qr/cutrank.json', 'w'), indent=0)
print('MAX DROP', max(r['drop'] for r in res))
print('PASS_EXPQR_CUTRANK' if res and max(r['drop'] for r in res) == 0 else 'FAIL_EXPQR_CUTRANK')   # F1 dead: drop 0 < 0.02 W on every cut
