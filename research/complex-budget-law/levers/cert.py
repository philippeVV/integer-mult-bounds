"""expBK obstruction certificate for cone repair of the 34 scout near-miss helpers (script-checkable, mod-p linear algebra).

Rewrite class certified impossible (per helper h, every hyperplane U = y^perp, hence every proper subspace):
  (i)  h's own A/B gates may be re-associated: each moved anywhere inside its exact syntactic-commutation window
       (h-gate order kept; earliest-feasible greedy is optimal since frames only grow);  U must contain every frame
       at which h is touched (own-frame condition);
  (ii) after h's last B gate, ANY number of compensating gates  h += c*g'  (g' not a target, not h) may be inserted at
       any B positions tau, each at frame  D = F_h(tau) + F_g'(tau) which must lie in U; the frame raise this forces
       on g' propagates forward through every later gate (monotone closure) and must respect every pin: targets'
       port-level hyperplanes, the terminal frames of already low-ended helpers, the frames of their early-cleanup tail
       gates.  Pin(g', tau) is computed exactly by a backward intersection DP; the condition is D <= Pin(g', tau);
  goal: row h of P^-1 restricted to OUT(U) = {g : end-of-B frame of g not inside U} must vanish.  Compensations change
       row h by -c * (row g' of S_tau^-1), S_tau = B-suffix map, so success needs  w|OUT in span{ row g' of S_tau^-1 |OUT }
       over the admissible (g', tau).  Raises are by vectors inside U, so OUT(U) does not change.
  Non-membership mod p (p odd, all coefficients dyadic) implies non-membership over Q, i.e. impossibility.
Also reports the 'own-frame' hyperplanes: for the scout's single-offender hyperplanes U, (i) is already infeasible.
usage: cert.py RAW.json.gz CENSUS.pkl [h ...]"""
import sys, json, pickle, bisect
import numpy as np
from collections import Counter
from common import *

PR = 2147483629
word = load(sys.argv[1]); d = pickle.load(open(sys.argv[2], 'rb'))
P, Pinv, endB0, touches, start = d['P'], d['Pinv'], d['endB'], d['touches'], d['start']
v = word['report']['v']; n = 2*v + word['report']['R']
B = word['B']
NOPINS = '--nopins' in sys.argv
hs = [int(x) for x in sys.argv[3:] if not x.startswith('--')] or sorted(int(h) for h in json.load(open('singles.json')))
FULLB = basis([1 << i for i in range(9)])


def ann(F):
    return basis([y for y in range(1, 512) if all(par(x & y) == 0 for x in F)])


def commute(a, b):
    ra, rb = set(regs(a)), set(regs(b))
    if a[1] not in rb and b[1] not in ra: return True
    return a[0] == 'add' and b[0] == 'add' and a[1] == b[1] and a[2] != b[1] and b[2] != a[1]


def inv(x):
    return pow(int(x) % PR, PR - 2, PR)


def fr(q):  # Fraction-like (num, den) -> mod p
    return q[0] % PR * inv(q[1]) % PR


# end-of-B pins (annihilator bases)
pin_end = {}
ONLYT, ONLYH = '--onlytargets' in sys.argv, '--onlyhpe' in sys.argv
for i in range(v):
    if not ONLYH: pin_end[v + i] = basis([word['ports'][i]])
for r, f in word.get('helper_final', {}).items():
    if not ONLYT: pin_end[int(r)] = ann(basis(f))
for i, f in ([] if ONLYT else word.get('early_cleanup', [])):
    e = word['cleanup'][i]
    for r in regs(e):
        pin_end[r] = basis(list(pin_end.get(r, ())) + list(ann(basis(f))))

Atouch = {r: [F for p, i, F in touches[r] if p == 'A'] for r in range(n)}


def placements(h, y):
    """greedy earliest placement of h's B gates inside y^perp (see ownfix.py); returns new B or None."""
    Btouch = {r: [(i, F) for p, i, F in touches[r] if p == 'B'] for r in range(n)}
    hB = [i for i, F in Btouch[h]]; hBset = set(hB)
    F0 = basis([x for F in Atouch[h] for x in F] + list(start[h]))
    if not in_yperp(F0, y): return None
    def fb(r, j):
        L = Btouch[r]; k = bisect.bisect_left([i for i, F in L], j)
        if k: return L[k - 1][1]
        return Atouch[r][-1] if Atouch[r] else start[r]
    Fh = F0; prev = -1; pos = []
    for i in hB:
        e = B[i]; lo = i
        while lo > 0 and (lo - 1) not in hBset and commute(B[lo - 1], e): lo -= 1
        hi = i
        while hi + 1 < len(B) and (hi + 1) not in hBset and commute(B[hi + 1], e): hi += 1
        other = [r for r in regs(e) if r != h]; found = None
        for j in range(max(lo, prev + 1), hi + 1):
            F = basis(list(Fh) + [x for r in other for x in fb(r, j)])
            if in_yperp(F, y): found = (j, F); break
            if not all(in_yperp(fb(r, j), y) for r in other): break
        if found is None: return None
        prev, Fh = found; pos.append((i, prev))
    # build B': move each h gate i to just before original position j (j <= i or j > i), then monotone closure
    moved = {i: j for i, j in pos}
    seq = []
    for k, e in enumerate(B):
        if k in moved: continue
        seq.append((k, 0, e))
    for i, j in moved.items():
        seq.append((j - 0.5 if j <= i else j + 0.5, 1, B[i]))
    seq.sort(key=lambda t: (t[0], t[1]))
    nb = [e for _, _, e in seq]
    cur = [None] * n
    for r in range(n):
        cur[r] = Atouch[r][-1] if Atouch[r] else start[r]
    out = []
    for e in nb:
        F = basis(list(e[-1]) + [x for r in regs(e) for x in cur[r]])
        for r in regs(e): cur[r] = F
        out.append(e[:-1] + [list(F)])
    return out, cur


def certify(h, y, nb, endB):
    U_ok = lambda F: in_yperp(F, y)
    OUT = [g for g in range(n) if endB[g] is not None and not U_ok(endB[g])]
    col = {g: k for k, g in enumerate(OUT)}
    m = len(OUT)
    w = np.zeros(m, dtype=np.int64)
    for g, x in Pinv[h].items():
        if g in col: w[col[g]] = x.numerator % PR * inv(x.denominator) % PR
    if not w.any():
        return dict(OUT=m, need=False)
    # backward sweep
    R = {}
    for g in OUT:
        z = np.zeros(m, dtype=np.int64); z[col[g]] = 1; R[g] = z
    pins = dict(pin_end)
    # frames before position k in nb
    pos = {}
    for k, e in enumerate(nb):
        for r in regs(e): pos.setdefault(r, []).append(k)
    def fb(r, k):
        L = pos.get(r, []); j = bisect.bisect_left(L, k)
        if j: return basis(nb[L[j - 1]][-1])
        return Atouch[r][-1] if Atouch[r] else start[r]
    lastH = max(pos[h])
    Fh_after = basis(nb[lastH][-1])
    basisV = {}   # pivot -> vector (mod p row echelon)
    def add(vec):
        vec = vec.copy()
        for piv in sorted(basisV):
            if vec[piv]:
                vec = (vec - vec[piv] * basisV[piv]) % PR
        nz = np.nonzero(vec)[0]
        if len(nz):
            p0 = nz[0]; vec = vec * inv(vec[p0]) % PR
            for q in list(basisV):
                if basisV[q][p0]: basisV[q] = (basisV[q] - basisV[q][p0] * vec) % PR
            basisV[p0] = vec
    def admissible(g, k):
        if g == h or v <= g < 2*v: return False
        Fg = fb(g, k)
        if not U_ok(Fg): return False
        D = basis(list(Fh_after) + list(Fg))
        return NOPINS or all(par(x & a) == 0 for x in D for a in pins.get(g, ()))
    nadm = 0
    cand = set(range(n))
    for k in range(len(nb), lastH, -1):
        # candidates at position k (insert before nb[k]); rows R are rows of S_k^-1 restricted to OUT
        for g in cand:
            if g in R and R[g].any() and admissible(g, k):
                add(R[g]); nadm += 1
        if k - 1 <= lastH: break
        e = nb[k - 1]
        # row update: S_{k-1}^-1 = M^-1 S_k^-1 ; gate t += c s  ->  row t -= c row s ; scale t by q -> row t /= q
        if e[0] == 'add':
            t, s = e[1], e[2]; c = e[3] % PR * inv(e[4]) % PR
            if s in R:
                R[t] = ((R[t] if t in R else np.zeros(m, dtype=np.int64)) - c * R[s]) % PR
        else:
            t = e[1]
            if t in R: R[t] = R[t] * (e[3] % PR) % PR * inv(e[2]) % PR
        # pin DP: both registers inherit the intersection
        a = basis(list(pins.get(e[1], ())) + (list(pins.get(e[2], ())) if e[0] == 'add' else []))
        for r in regs(e): pins[r] = a
        cand = set(regs(e))
    # membership test
    vec = w.copy()
    for piv in sorted(basisV):
        if vec[piv]: vec = (vec - vec[piv] * basisV[piv]) % PR
    return dict(OUT=m, need=True, dimV=len(basisV), admissible_vectors=nadm, w_in_span=not vec.any())


res = {}
for h in hs:
    rows = []
    feas_y = 0
    for y in range(1, 512):
        pl = placements(h, y)
        if pl is None: continue
        nb, curend = pl
        feas_y += 1
        endB = list(curend)
        r = certify(h, y, nb, endB)
        r['y'] = y
        rows.append(r)
    ok = [r for r in rows if (not r['need']) or r['w_in_span']]
    res[h] = dict(own_feasible_hyperplanes=feas_y, repairable_hyperplanes=[r['y'] for r in ok],
                  max_dimV=max((r.get('dimV', 0) for r in rows), default=0))
    print(h, json.dumps(res[h]), flush=True)
json.dump({str(k): v for k, v in res.items()}, open('cert' + ''.join(a for a in sys.argv if a.startswith('--')) + '.json', 'w'), indent=0)
print('helpers with a repairable hyperplane in this class:', sum(1 for r in res.values() if r['repairable_hyperplanes']), 'of', len(res))
