"""(B'') Re-ordering-invariant cone bound: any re-ordering of the gate list that keeps non-commuting gates in
order (Mazurkiewicz trace of forced.build, compensation kept as real gates) AND any frame refit.

For framed gates p, g sharing a register with p < g in the dependency order, every schedule runs p before g,
so F_p <= F_g. Hence
    lo(g) = span(source port lines of g's sources (outside warps), helper/target start pins, lo(p) for such p)
    hi(g) = meet(target hyperplanes / end pins of g's registers, hi(s) for every such successor s)
are valid frame bounds for all schedules and refits (no A/B boundary pins are used: they are order-specific).
Order-free stream bound: in any schedule the first frame of stream h is <= hi(g) and the last >= lo(g) for
all its gates g, so the touch-span rank is >= max_g dim lo(g) - min_g dim hi(g) =: w_h and (subadditivity
of f(r) = r ln(73/r)) its touch-span log-mass is >= f(w_h).
usage: propdag.py WORD.json.gz
"""
import sys, json, math
from collections import Counter, defaultdict
from wordmodel import Word, basis, contained, regs
from propbound import meet, span, f
import forced


def main():
    Wd = Word(sys.argv[1])
    w = Wd.w; v = Wd.v; FULL = Wd.full; ports = w['ports']
    keep, idx, N, desc, anc, streams = forced.build(Wd, 'literal', False)
    ev = Wd.events
    warped_at, wts = {}, defaultdict(list)
    for t, (s, e) in enumerate(ev):
        if e[0] == 'warp':
            warped_at.setdefault(e[1], t)
            wts[e[1]].append(t)
    framed = [t for t in keep if t in Wd.gateframe]
    onreg = defaultdict(list)
    for t in framed:
        for q in regs(ev[t][1]):
            onreg[q].append(t)

    def slo(q, t):
        if q in warped_at:
            return ()
        if q < v:
            return basis([ports[q]])
        return () if v <= q < 2 * v else Wd.start[q]

    def shi(q, t):
        if q in warped_at:
            return FULL
        return Wd.end[q]
    # warp episodes: no frame relation across a warp on the lender
    def same_epoch(q, a, b):
        return q not in warped_at or sum(1 for x in wts[q] if x < a) == sum(1 for x in wts[q] if x < b)
    lo, hi = {}, {}
    for t in framed:
        i = idx[t]
        L = []
        for q in regs(ev[t][1]):
            L.extend(slo(q, t))
            for p in onreg[q]:
                if p < t and (anc[i] >> idx[p]) & 1 and same_epoch(q, p, t):
                    L.extend(lo[p])
        lo[t] = basis(L)
    for t in reversed(framed):
        i = idx[t]
        H = FULL
        for q in regs(ev[t][1]):
            H = meet(H, shi(q, t))
            for s in onreg[q]:
                if s > t and (desc[i] >> idx[s]) & 1 and same_epoch(q, s, t):
                    H = meet(H, hi[s])
        hi[t] = H
    ctrl = sorted(t for t in framed if not (contained(lo[t], Wd.gateframe[t]) and contained(Wd.gateframe[t], hi[t])))
    ctrl = sorted(ctrl, key=lambda t: -t)
    if ctrl:
        t = ctrl[0]; print("CTRL", t, ev[t], lo[t], hi[t], Wd.gateframe[t], file=sys.stderr)
        for q in regs(ev[t][1]):
            print("  reg", q, [(s, ev[s][1], hi[s]) for s in onreg[q] if s > t and hi[s] == ()][:3], file=sys.stderr)
        bad0 = min(ctrl); print("first bad", bad0, ev[bad0], lo[bad0], hi[bad0], Wd.gateframe[bad0], file=sys.stderr)
    assert not ctrl
    rows = []
    for kind, r, m in streams:
        if kind == 'source':
            continue
        ts = [keep[i] for i in range(N) if (m >> i) & 1 and keep[i] in Wd.gateframe]
        if not ts:
            continue
        wv = max(len(lo[t]) for t in ts) - min(len(hi[t]) for t in ts)
        dims = []
        for t in ts:
            d = len(Wd.gateframe[t])
            if not dims or d != dims[-1]:
                dims.append(d)
        act = sum(f(b - a) for a, b in zip(dims, dims[1:]))
        rows.append(dict(reg=r, kind=kind, w=max(wv, 0), lb=f(max(wv, 0)), actual=act, actual_rank=dims[-1] - dims[0]))
    out = dict(word=sys.argv[1], streams=len(rows), span_actual=sum(x['actual'] for x in rows),
               lb_orderfree=sum(x['lb'] for x in rows),
               rank_actual=sum(x['actual_rank'] for x in rows), rank_lb=sum(x['w'] for x in rows),
               w_hist=dict(sorted(Counter(x['w'] for x in rows).items())),
               gate_slack_hist=dict(sorted(Counter(len(hi[t]) - len(lo[t]) for t in framed).items())))
    out['ratio'] = out['lb_orderfree'] / out['span_actual']
    # price the envelope: every helper register keeps only one climb of rank w_h (entry/exit free)
    prof = Counter(Wd.paid)
    for r in Wd.helpers:
        for c in Wd.climbs.get(r, []):
            prof[c[1]] -= 1
    for x in rows:
        if x['w'] > 0:
            prof[x['w']] += 1
    prof = +prof
    sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', 'complex-local-word-ceiling'))
    import ceiling as CE
    from fractions import Fraction as Q
    g = Q(1, 10 ** 12)
    b0 = CE.ceiling(dict(Wd.paid), 73, 896, 9, 64, grid=g)[1]
    b1 = CE.ceiling(dict(prof), 73, 896, 9, 64, grid=g)[1]
    out['envelope'] = dict(profile=dict(sorted(prof.items())), ceiling_ratio=float(b1 / b0))
    if w['warps']:
        from chains import price_root
        out['envelope']['t74_root'] = price_root(dict(Wd.paid), dict(prof))
    # all helper climbs removed (frame-local ideal) for comparison
    prof0 = Counter(Wd.paid)
    for r in Wd.helpers:
        for c in Wd.climbs.get(r, []):
            prof0[c[1]] -= 1
    out['ideal_ceiling_ratio'] = float(CE.ceiling(dict(+prof0), 73, 896, 9, 64, grid=g)[1] / b0)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
