import os
"""Frame model of a BORROW-WORD replicating audit_local.run_borrow's frame chain (section 4), with minimal forward
lifting from raised start frames.  Events: chronological list of (kind, regs, frame, ref) where ref locates the item."""
import json, gzip, heapq
from collections import Counter, defaultdict

def basis(vs):
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

def contained(a, b): return basis(list(a) + list(b)) == b

def gmul(cols, x):
    y = 0; i = 0
    while x:
        if x & 1: y ^= cols[i]
        x >>= 1; i += 1
    return y

def registers(e): return e[1:3] if e[0] == 'add' else e[1:2]

class Word:
    def __init__(self, w):
        self.w = w
        v, h, R = w['report']['v'], w['report']['h'], w['report']['R']
        self.v, self.h, self.n = v, h, 2 * v + R
        n = self.n
        full = basis([1 << i for i in range(h)]); self.full = full
        self.removed = set(w['removed'])
        self.live = [r for r in range(n) if r not in self.removed]
        end = [full] * v + [basis([x for x in range(1, 1 << h) if (x & p).bit_count() % 2 == 0]) for p in w['ports']] + [full] * R
        for q, U in w['helper_final'].items(): end[int(q)] = basis(U)
        start = [(p,) for p in w['ports']] + [()] * (v + R)
        for q, U in w['start_extra'].items(): start[int(q)] = basis(U)
        self.end, self.start = end, start
        ev = []   # (kind, regs, frame, ref)
        npre = w.get('prelude_events', 0)
        for i, e in enumerate(w['A'][:npre]):
            if e[0] == 'warp':
                ev.append(('pwarp', (e[1],), (basis(e[3]), basis(e[4])), ('A', i)))
            else:
                ev.append(('gate', tuple(registers(e)), basis(e[-1]), ('A', i)))
        for i in range(npre, len(w['A'])):
            e = w['A'][i]; ev.append(('gate', tuple(registers(e)), basis(e[-1]), ('A', i)))
        ev.append(('checkA', (), None, None))
        bypos = defaultdict(list)
        for op, rec in enumerate(w['warps'], 1): bypos[rec[0]].append(op)
        B = w['B']
        for p in range(len(B) + 1):
            for op in bypos.get(p, []):
                rec = w['warps'][op - 1]
                ev.append(('warp', (rec[1],), (basis(rec[3]), basis(rec[4])), ('W', op - 1)))
            if p < len(B):
                e = B[p]; ev.append(('gate', tuple(registers(e)), basis(e[-1]), ('B', p)))
        for i, e in enumerate(w['restore']): ev.append(('gate', tuple(registers(e)), basis(e[-1]), ('restore', i)))
        self.ev = ev
        self.byreg = defaultdict(list)
        for k, (kind, regs, fr, ref) in enumerate(ev):
            for r in regs: self.byreg[r].append(k)
        self.lev = {int(t): basis(U) for t, U in w.get('target_levels', {}).items()}

    def chain(self, r, frames=None, start=None):
        """list of frames of register r: start, its events..., end (warps: two entries, the warp resets)"""
        frames = frames or {}
        out = [[start if start is not None else self.start[r]]]
        for k in self.byreg[r]:
            kind, regs, fr, ref = self.ev[k]
            f = frames.get(k, fr)
            if kind in ('warp', 'pwarp'):
                out[-1].append(f[0]); out.append([f[1]])
            else: out[-1].append(f)
        out[-1].append(self.end[r])
        return out

    def climbs(self, r, frames=None, start=None):
        c = []
        for seg in self.chain(r, frames, start):
            for a, b in zip(seg, seg[1:]):
                d = len(b) - len(a)
                if d: c.append(d)
        return c

    def lift(self, newstart):
        """newstart: {reg: basis}.  Minimal forward lift.  Returns (frames dict ev->new frame, flags, touched regs)."""
        cur = {}       # reg -> lifted current frame (only regs that differ from original)
        frames = {}; flags = set(); touched = set(newstart)
        heap = []
        pos = {}
        for r, U in newstart.items():
            cur[r] = U
            ks = self.byreg[r]
            if ks: heapq.heappush(heap, ks[0])
            else: flags.add('noevents')
        seen = set()
        # original frame of register r just before event k is needed when partner not lifted: that's ev frame of its prev event; we
        # only need lifted ones (non-lifted partners are at their original frame which is <= gate frame).
        while heap:
            k = heapq.heappop(heap)
            if k in seen: continue
            seen.add(k)
            kind, regs, fr, ref = self.ev[k]
            if kind == 'gate':
                U = fr
                extra = [cur[r] for r in regs if r in cur]
                U2 = basis(list(U) + [x for b in extra for x in b])
                if U2 != U:
                    frames[k] = U2
                    if ref[0] == 'A' and ref[1] < self.w.get('prelude_events', 0): flags.add('prelude-gate')
                for r in regs:
                    touched.add(r)
                    if U2 != U:
                        cur[r] = U2
                    else:
                        cur.pop(r, None)
                    nxt = self._next(r, k)
                    if r in cur and nxt is not None: heapq.heappush(heap, nxt)
                    if r in cur and nxt is None:
                        if not contained(cur[r], self.end[r]): flags.add('end:%d' % r)
                        if self.v <= r < 2 * self.v: flags.add('target-end')
            elif kind in ('warp', 'pwarp'):
                r = regs[0]; U, E = fr
                if r in cur and cur[r] != U:
                    rec = self.w['warps'][ref[1]] if kind == 'warp' else self.w['warps'][abs(self.w['A'][ref[1]][2]) - 1]
                    cols = rec[2]
                    if kind == 'pwarp' and self.w['A'][ref[1]][2] < 0: flags.add('pwarp-inverse')
                    E2 = basis([gmul(cols, x) for x in cur[r]])
                    frames[k] = (cur[r], E2); flags.add('warp')
                    cur[r] = E2
                    nxt = self._next(r, k)
                    if nxt is not None: heapq.heappush(heap, nxt)
                    else:
                        if not contained(cur[r], self.end[r]): flags.add('end:%d' % r)
                else: cur.pop(r, None)
            elif kind == 'checkA':
                pass
        # registers lifted whose next event is the end and never processed: newstart regs w/o events handled.
        # A-end checks: targets and centers must be at their original frame at the end of A (index of checkA)
        ka = next(i for i, e in enumerate(self.ev) if e[0] == 'checkA')
        for r in touched:
            if self.v <= r < 2 * self.v or r in self.w['centers']:
                # frame at end of A: last event of r before ka
                prev = [k for k in self.byreg[r] if k < ka]
                if prev and prev[-1] in frames: flags.add('A-end:%d' % r)
                if not prev and r in newstart: flags.add('A-end:%d' % r)
        return frames, flags, touched

    def _next(self, r, k):
        ks = self.byreg[r]
        import bisect
        i = bisect.bisect_right(ks, k)
        return ks[i] if i < len(ks) else None

    def profile(self):
        p = Counter()
        for r in self.live:
            for d in self.climbs(r): p[d] += 1
        return p

def load(path):
    return json.loads(gzip.decompress(open(path, 'rb').read()))
