"""expNV / stage.py -- can a second centre call reuse a fresh stage in the retained rho_(s,q) routes?
Exact GF(2) linear algebra on the actual compressed window V = F2^9 + (F2^8)^8 (m = 73), the 112 archived
port labels (DATA-PHASES.json, same file the fivecall audit pins) and the inherited U_q / rho_(s,q) formula
(INHERITED-COMPLEX-GEOMETRY.txt).  Stdlib only.  Run: python3 -I -B stage.py

Admission condition of a centre call on a data family with binary frame W (the span on which C acts):
the routed local view must be the call's input frame <q>, i.e.  rho(W) cap L9 = <q>  (L9 = local nine-space),
and the call then adds q^perp locally, i.e. one exterior 8-block globally (nested transition, rank 8).
Checks:
 A  fresh stage: W_{s-1} = L9 + E_1..E_{s-1}: rho_(s,q)W cap L9 = <q> for all 112 ports, s = 1..8 (admitted),
    and the call output pulled back is W_s (rank +8).
 B  shared stage (same family, stage s' <= current level s, ANY port p, and any p->q->p excursion
    rho_(s',p) rho_(s',q)): the routed local view is all of L9 (dimension 9), so no call is admissible.
 C  route independence: for ANY g in O(73,2), dim(gW cap L9) >= dim W - 64, so a frame of dimension > 65
    never admits a call; the cheapest repair is a nested descent of Fourier rank dim W - 65 (= 8 at W = F2^73),
    i.e. one full centre rank.  Checked against random orthogonal-ish products of the actual rho's as a control.
"""
import json, random, sys
import os
from pathlib import Path
ROOT = Path(os.environ.get('IMB_REPO') or Path(__file__).resolve().parents[4])
ROUTE = str(ROOT / 'research/direct-delivery504-global/inputs/data-base/work/complex-nine-native/DATA-PHASES.json')
M = 73; L9 = (1 << 9) - 1
def blk(s): return ((1 << 8) - 1) << (9 + 8 * (s - 1))
def dot(x, y): return (x & y).bit_count() & 1
def u_columns(q):
    j = next(k for k in range(9) if not (q >> k) & 1); v = q ^ (1 << j); cols = []
    for i in range(9):
        x = 1 << (j if i == 0 else 0 if i == j else i)
        if dot(x, v): x ^= v
        cols.append(x)
    assert cols[0] == q and all(dot(cols[i], cols[k]) == (i == k) for i in range(9) for k in range(9))
    return cols
COL = {}
def rho(q, x, s):
    """rho_(s,q): fixes q and all other exterior blocks, exchanges q^perp (cols 1..8 of U_q) with E_s."""
    c = COL[q]; sh = 9 + 8 * (s - 1)
    loc = x & L9; fr = (x >> sh) & 255; out = x & ~(L9 | (255 << sh))
    nl = q if dot(loc, q) else 0
    for i in range(8):
        if (fr >> i) & 1: nl ^= c[i + 1]
    nf = sum(dot(loc, c[i + 1]) << i for i in range(8))
    return out | nl | (nf << sh)
def basis(vs):
    b = {}
    for x in vs:
        for p in sorted(b, reverse=True):
            if x >> p & 1: x ^= b[p]
        if x:
            p = x.bit_length() - 1
            for k in b:
                if b[k] >> p & 1: b[k] ^= x
            b[p] = x
    return sorted(b.values())
def span_bits(mask): return [1 << i for i in range(M) if mask >> i & 1]
def inter_dim(A, B): return len(basis(A)) + len(basis(B)) - len(basis(A + B))
def local_view(W):
    """basis of W cap L9 (W given by a basis)."""
    # W cap L9 = vectors of W with zero exterior part: kernel of projection to exterior
    ext = [w >> 9 for w in W]; rows = list(zip(W, ext)); piv = {}
    kernel = []
    for w, e in rows:
        comb = w
        for p in sorted(piv, reverse=True):
            if e >> p & 1: e ^= piv[p][1]; comb ^= piv[p][0]
        if e: piv[e.bit_length() - 1] = (comb, e)
        else: kernel.append(comb)
    return basis(kernel)
def main():
    ports = sorted({r['q'] for r in json.load(open(ROUTE))['recipes']})
    assert len(ports) == 112 and all(q.bit_count() % 2 for q in ports)
    for q in ports: COL[q] = u_columns(q)
    # involution / orthogonality spot checks on the 73-bit window
    for q in ports[:8]:
        for s in (1, 8):
            assert all(rho(q, rho(q, 1 << i, s), s) == 1 << i for i in range(M))
    Wlev = [L9 | sum(blk(t) for t in range(1, s + 1)) for s in range(9)]   # W_s, s = 0..8 (dim 9+8s)
    # A: fresh stage admits the call for every port and stage; output is W_s
    nA = 0
    for s in range(1, 9):
        W = span_bits(Wlev[s - 1])
        for q in ports:
            R = [rho(q, w, s) for w in W]
            lv = local_view(R)
            assert len(lv) == 1 and lv[0] == q, (s, q, lv)
            out = basis(R + span_bits(L9))                     # call: local frame <q> -> L9 (exterior factor kept)
            back = basis([rho(q, w, s) for w in out])
            assert len(back) == len(W) + 8 and basis(back + span_bits(Wlev[s])) == basis(span_bits(Wlev[s]))
            nA += 1
    print(f'A fresh stage: {nA} (stage,port) admissions, local view = <q>, output frame W_s (rank +8)  OK')
    # B: shared stage, same family: every s' <= s, every port p, and every p->q excursion
    nB = nX = 0
    for s in range(1, 9):
        W = span_bits(Wlev[s])
        for s2 in range(1, s + 1):
            for p in ports:
                R = [rho(p, w, s2) for w in W]
                assert len(local_view(R)) == 9; nB += 1
            for p, q in random.Random(s * 100 + s2).sample([(a, b) for a in ports for b in ports if a != b], 40):
                R = [rho(p, rho(q, w, s2), s2) for w in W]
                assert len(local_view(R)) == 9; nX += 1
    print(f'B shared stage: {nB} (level,stage,port) routes and {nX} p->q excursions all give local view = L9 '
          f'(dim 9 != 1): no call admissible without repair  OK')
    # C: route independence and repair rank
    rng = random.Random(7); full = span_bits((1 << M) - 1)
    for trial in range(200):
        W = span_bits(Wlev[rng.randrange(8, 9)])                 # dim 73
        g = [(rng.choice(ports), rng.randrange(1, 9)) for _ in range(rng.randrange(1, 12))]
        R = W
        for q, s in g: R = [rho(q, w, s) for w in R]
        assert len(local_view(R)) >= len(W) - 64
    # minimal descent: a subframe W' of W_8 (codim k) admits a call iff dim W' <= 65 ... and W' = W_8 minus one block works
    for s in range(1, 9):
        Wp = span_bits(Wlev[8] & ~blk(s)); assert len(Wp) == 65
        q = ports[s]; R = [rho(q, w, s) for w in Wp]
        assert local_view(R) == [q]
    print('C any route g: dim(gW cap L9) >= dim W - 64 (200 random rho-words on W_8); at W_8 (dim 73) the cheapest '
          're-admitting repair is a nested descent of rank 8 (drop one block), which then admits the call at that stage  OK')
    print('VERDICT: a second call on the same family at a used stage needs a paid repair of rank >= 8 per port family '
          '(= one full centre rank): stop test FIRES for same-family stage sharing')
    print('PASS_EXPNV_STAGE')
if __name__ == '__main__':
    main()
