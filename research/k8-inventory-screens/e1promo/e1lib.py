"""E1 -> E2 entrance-promotion analysis library (expB).  Exact GF(2) subspace arithmetic via 512-bit membership masks.

Mirrors research/gauge44-leaver8/independent/audit_local.py's frame/chain model exactly:
  start: sources (port,), targets (), helpers () or start_extra
  climbs: A gates, [boundary: targets at target_levels, centres at centre frames], B gates, [endB],
          early-cleanup tail gates (frame from early_cleanup), end climbs (sources full, targets port hyperplane,
          helpers full or helper_final)
  cones: for each low-end helper q, endB[k] must lie in end[q] for k in support(pinv row q) and dirt column q.
All legality conditions are upper bounds on frames, so the minimal monotone completion of raised start frames is legal iff
each raised start lies in its backward ceiling (computed here exactly).
"""
import gzip, json, hashlib, math, functools
from fractions import Fraction as F
from collections import Counter

FULLMASK = (1 << 512) - 1


def regs(e):
    return list(e[1:3]) if e[0] == 'add' else list(e[1:2])


def basis(vs):  # identical canonical RREF to audit_local.basis
    piv = {}
    for a in vs:
        while a:
            bit = a.bit_length() - 1
            if bit in piv: a ^= piv[bit]
            else:
                piv[bit] = a; break
    for k in sorted(piv):
        for j in piv:
            if j > k and piv[j] >> k & 1: piv[j] ^= piv[k]
    return tuple(piv[k] for k in sorted(piv, reverse=True))


@functools.lru_cache(maxsize=None)
def mask_of_basis(b):
    m = 1
    for x in b:
        if not (m >> x) & 1:
            new = m
            mm = m
            while mm:
                low = mm & -mm; y = low.bit_length() - 1; mm ^= low
                new |= 1 << (y ^ x)
            m = new
    return m


def mask(vs):
    return mask_of_basis(basis(list(vs)))


def members(m):
    out = []
    while m:
        low = m & -m; out.append(low.bit_length() - 1); m ^= low
    return out


@functools.lru_cache(maxsize=None)
def basis_of_mask(m):
    return basis(members(m))


def dim(m):
    return bin(m).count('1').bit_length() - 1


@functools.lru_cache(maxsize=None)
def add(m1, m2):
    if m1 & ~m2 == 0: return m2
    if m2 & ~m1 == 0: return m1
    return mask_of_basis(basis(list(basis_of_mask(m1)) + list(basis_of_mask(m2))))


def load(path):
    raw = gzip.decompress(open(path, 'rb').read())
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


class Word:
    def __init__(self, w):
        self.w = w
        v, R = w['report']['v'], w['report']['R']; self.v, self.R = v, R; n = self.n = 2 * v + R
        h = w['report']['h']; assert h == 9
        self.full = mask(range(1 << 9)) if False else mask([1 << i for i in range(9)])
        ev = []
        for i, e in enumerate(w['A']): ev.append(('A', i, regs(e), mask(e[-1])))
        self.nA = len(ev)
        for i, e in enumerate(w['B']): ev.append(('B', i, regs(e), mask(e[-1])))
        self.nB = len(ev)
        for j, (i, U) in enumerate(w.get('early_cleanup', [])): ev.append(('T', j, regs(w['cleanup'][i]), mask(U)))
        self.ev = ev
        self.start = [mask([p]) for p in w['ports']] + [1] * (v + R)
        for q, U in w.get('start_extra', {}).items(): self.start[int(q)] = mask(U)
        self.end = [self.full] * v + [mask([x for x in range(1, 512) if (x & p).bit_count() % 2 == 0]) for p in w['ports']] + [self.full] * R
        for q, U in w.get('helper_final', {}).items(): self.end[int(q)] = mask(U)
        self.entering = sorted(int(q) for q in w.get('start_extra', {}))
        self.lowends = sorted(int(q) for q in w.get('helper_final', {}))
        # boundary ceilings (exact frames at A/B boundary)
        lev = {int(t): mask(U) for t, U in w.get('target_levels', {}).items()}
        self.boundary = {v + t: lev.get(t, 1) for t in range(v)}
        for q, U in zip(w['centers'], w['center_frames']): self.boundary[q] = mask(U)
        # register event lists
        self.revs = [[] for _ in range(n)]
        for gi, (_, _, rs, _) in enumerate(ev):
            for q in rs: self.revs[q].append(gi)
        self.cone = None

    def forward(self, start=None, frames=None):
        """frames before each event per register; returns (cur_before[(gi,q)], endB, final cur, paid Counter excl. centre paths)."""
        start = self.start if start is None else start
        fr = [e[3] for e in self.ev] if frames is None else frames
        cur = list(start); paid = Counter(); before = {}
        endB = None
        for gi, (sec, i, rs, _) in enumerate(self.ev):
            if gi == self.nA:
                for q, U in self.boundary.items(): assert cur[q] == U, ('boundary', q)
            if gi == self.nB: endB = list(cur)
            U = fr[gi]
            for q in rs:
                assert cur[q] & ~U == 0, ('nonmonotone', gi, q)
                before[(gi, q)] = cur[q]
                d = dim(U) - dim(cur[q])
                if d: paid[d] += 1
                cur[q] = U
        if endB is None: endB = list(cur)
        for q in range(self.n):
            assert cur[q] & ~self.end[q] == 0, ('end', q)
            d = dim(self.end[q]) - dim(cur[q])
            if d: paid[d] += 1
        return before, endB, cur, paid

    def compute_cones(self, cache=None):
        """support sets for low-end cones, exact rationals (mirrors audit_local.run)."""
        if cache:
            try:
                d = json.load(open(cache))
                if d['sha'] == self.sha:
                    self.cone = {int(k): v for k, v in d['cone'].items()}; return self.cone
            except Exception:
                pass
        w = self.w; v, n = self.v, self.n

        def ident(): return [{k: F(1)} for k in range(n)]

        def ap(e, rows):
            if e[0] == 'add':
                t, s = e[1:3]; q = F(e[3], e[4]); dst = rows[t]
                for k, x in list(rows[s].items()):
                    a = dst.get(k, F(0)) + q * x
                    if a: dst[k] = a
                    else: dst.pop(k, None)
            else:
                t = e[1]; q = F(e[2], e[3]); rows[t] = {k: q * x for k, x in rows[t].items()}

        def inv(e):
            r = list(e)
            if r[0] == 'add': r[3] = -r[3]
            else: r[2], r[3] = r[3], r[2]
            return r
        main = w['A'] + w['cut'] + w['B']
        raw = ident()
        for e in main: ap(e, raw)
        dirt = [{} for _ in range(n)]
        for y in range(v, 2 * v): dirt[y] = {k: a for k, a in raw[y].items() if k >= 2 * v}
        for e in w['precondition']: ap(e, dirt)
        comp = [['add', y, k, -a.numerator, a.denominator, []] for y in range(v, 2 * v) for k, a in sorted(dirt[y].items())]
        self.compensation_targets = sorted({e[2] for e in comp})
        prefix = w['precondition'] + comp + main
        p = ident(); pinv = ident()
        for e in prefix: ap(e, p)
        for e in reversed(prefix): ap(inv(e), pinv)
        cone = {}
        for q in self.lowends:
            K = set(pinv[q]) | {k for k in range(n) if k != q and p[k].get(q)}
            cone[q] = sorted(K)
        self.cone = cone
        if cache: json.dump(dict(sha=self.sha, cone={str(k): v for k, v in cone.items()}, comp_cols=self.compensation_targets), open(cache, 'w'))
        return cone

    def ceilings(self):
        """backward ceilings: Cstart[q] = largest frame q's start may have; gate_ceil[gi]."""
        n = self.n
        C = list(self.end)
        coneceil = [FULLMASK] * n
        for q, K in self.cone.items():
            for k in K: coneceil[k] &= self.end[q]
        self.coneceil = coneceil
        gc = [None] * len(self.ev)
        for gi in range(len(self.ev) - 1, -1, -1):
            sec, i, rs, U = self.ev[gi]
            # entering the B/T boundary going backward: apply cone ceilings once
            if gi == self.nB - 1:
                for k in range(n): C[k] &= coneceil[k]
            if gi == self.nA - 1:
                for q, U0 in self.boundary.items(): C[q] &= U0
            c = FULLMASK
            for q in rs: c &= C[q]
            gc[gi] = c
            for q in rs: C[q] = c
        # registers with no events at all: still subject to those constraints (handled by the order above only if events
        # exist on both sides; apply them here for completeness)
        self.gate_ceil = gc
        self.Cstart = C
        return C, gc

    def complete(self, newstart):
        """minimal monotone completion for raised starts {q: mask}; returns (raised gate frames {gi: mask}, dprofile Counter,
        legal flag, reasons, touched registers)."""
        cur = {}  # affected registers -> current new frame
        for q, m in newstart.items():
            assert self.start[q] & ~m == 0
            if m != self.start[q]: cur[q] = m
        raised = {}; dprof = Counter(); bad = []
        import heapq
        pos = {q: 0 for q in range(self.n)}  # pointer into revs (affected only)
        heap = []
        nxt = {}
        for q in cur:
            ev = self.revs[q]
            if ev: heapq.heappush(heap, ev[0])
        orig_before = self.orig_before
        processed = set()
        endB_new = {}
        state_at_boundary_checked = False

        def check_boundary(q, m):
            if q in self.boundary and m != self.boundary[q]: bad.append(('boundary', q))

        # registers whose first event is after the boundary / endB need checks when they pass those points
        def passes(q, gprev, gnext, m):
            # q holds frame m from after event gprev (index or -1) until before event gnext (index or len(ev))
            if gprev < self.nA <= gnext: check_boundary(q, m)
            if gprev < self.nB <= gnext:
                endB_new[q] = m
                if m & ~self.coneceil[q]: bad.append(('cone', q))
        for q, m in cur.items():
            ev = self.revs[q]
            passes(q, -1, ev[0] if ev else len(self.ev), m)
        touched = set(cur)
        while heap:
            gi = heapq.heappop(heap)
            if gi in processed: continue
            processed.add(gi)
            sec, i, rs, U = self.ev[gi]
            aff = [q for q in rs if q in cur]
            if not aff: continue
            Un = U
            for q in aff: Un = add(Un, cur[q])
            for q in rs:
                ob = orig_before[(gi, q)]
                nb = cur.get(q, ob)
                dold = dim(U) - dim(ob); dnew = dim(Un) - dim(nb)
                if dold != dnew:
                    if dold: dprof[dold] -= 1
                    if dnew: dprof[dnew] += 1
            if Un != U:
                raised[gi] = Un
                if Un & ~self.gate_ceil[gi]: bad.append(('gate_ceiling', gi))
                for q in rs:
                    cur[q] = Un; touched.add(q)
            else:
                for q in aff: del cur[q]
            for q in rs:
                ev = self.revs[q]; k = ev.index(gi)
                g2 = ev[k + 1] if k + 1 < len(ev) else len(self.ev)
                if q in cur:
                    passes(q, gi, g2, cur[q])
                    if g2 < len(self.ev): heapq.heappush(heap, g2)
        # end climbs
        for q, m in cur.items():
            if m & ~self.end[q]: bad.append(('end', q)); continue
            ev = self.revs[q]
            ob = self.final[q]
            dold = dim(self.end[q]) - dim(ob); dnew = dim(self.end[q]) - dim(m)
            if dold != dnew:
                if dold: dprof[dold] -= 1
                if dnew: dprof[dnew] += 1
        dprof = Counter({k: c for k, c in dprof.items() if c})
        return raised, dprof, not bad, bad, touched

    def centre_paid(self):
        w = self.w; v = self.v
        lev = {int(t): basis(U) for t, U in w.get('target_levels', {}).items()}
        paid = Counter()
        for q, U in zip(w['centers'], w['center_frames']):
            U = basis(U)
            served = {lev.get(e[1] - v, ()) for e in w['cut'] if e[2] == q}
            for T in sorted(served, key=lambda x: -len(x)):
                assert mask_of_basis(T) & ~mask_of_basis(U) == 0
                if T != U: paid[len(U) - len(T)] += 1
                U = T
        return paid

    def prepare(self):
        before, endB, final, paid = self.forward()
        self.orig_before = before; self.final = final; self.endB = endB; self.paid0 = paid
        return paid


# ---- pricing (caller-supplied frontier dual, approximate) ----
EXC = {1: 0.0072326, 2: 0.0121212, 3: 0.0161270, 4: 0.0195598, 5: 0.0225670, 8: 0.0297655}
BETA = math.log(1.0072326) / math.log(73)
for r in (6, 7, 9):
    EXC[r] = r * ((73 / r) ** BETA - 1)   # extrapolated from the same curve (not caller-supplied)
WPRICE = {4: -1, 5: -20, 7: -8, 8: -2, 9: 4}


def price(dprof, dwidths):
    """dprof: Counter rank->delta count; dwidths: Counter width->delta count. Returns (dRC, parts)."""
    dex = sum(EXC[r] * c for r, c in dprof.items())
    dM = sum(r * c for r, c in dprof.items())
    missing = [wd for wd in dwidths if dwidths[wd] and wd not in WPRICE]
    dP = sum(WPRICE.get(wd, 0) * c for wd, c in dwidths.items())
    return dex - 0.000784 * dM + 0.000293 * dP, dict(dExcess=dex, dM=dM, dP=dP, unpriced_widths=missing)
