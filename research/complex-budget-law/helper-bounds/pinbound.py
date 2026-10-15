"""(B) Static-pin lower bound on helper frame-chain log-mass (fixed gate list, any frame refit).

Frame calculus (audit_local): every framed gate runs at one frame F; each touched register climbs to F
(F contains its current frame) and frames only grow; register q always satisfies start(q) <= frame(q) <=
end(q). Static pins of a gate on helper h, from its partner q (valid for EVERY frame assignment of the
same gate list):
   lower:  start(q)  (+ <p_s> for a source s outside any warp episode)
   upper:  end(q)    (target t: the hyperplane p_t^perp; restored cold helper: helper_final; else FULL)
   centre registers additionally carry F = centre_frame at the cut (lower = upper = U_c).
Helper-helper couplings beyond the partners' static start/end pins are dropped (relaxation).

Lemma (static-pin chain bound). Let the framed gates touching stream h be g_1..g_k in time order with
lower pins L_i and upper pins U_i, start pin S, end pin E. In every frame refit of the fixed gate list,
   F_i >= S + L_1 + ... + L_i =: lo_i,   F_i <= E ^ U_i ^ U_{i+1} ^ ... ^ U_k =: hi_i,
so dim F_i lies in [dim lo_i, dim hi_i], nondecreasing, from dim S to dim E. Since f(r) = r ln(m/r)
is concave with f(0)=0, the minimum chain log-mass over such dimension sequences (an exact DP) is a lower
bound on h's frame-chain log-mass (and its minimal number of increments a lower bound on its climbs).
Proof: the true chain's dimension sequence is one of the sequences the DP ranges over. QED.

The script checks feasibility lo_i <= hi_i (as subspaces) on the actual word as a positive control,
compares LB with the actual per-stream log-mass (stream priced as a fresh register: start pin -> gate
frames -> end pin), and reports the 8D-read sub-lemma counts:
  (8D) if h is read into target t at a frame of dim 8 then that frame is p_t^perp, every later frame of h
       is p_t^perp or FULL, so h pays a rank-1 climb (8 -> 9, log-cost ln 73) unless E = p_t^perp, and no
       other target t' != t with p_t' != p_t can ever read h at or after that gate.

usage: pinbound.py WORD.json.gz [--json out]
"""
import sys, json, math
from collections import Counter, defaultdict
from wordmodel import Word, basis, contained, regs

M = 73


def f(r):
    return r * math.log(M / r) if r > 0 else 0.0


import functools
@functools.lru_cache(maxsize=None)
def meet(a, b, h=9):
    """intersection of subspaces given by bases (via annihilators)."""
    # annihilator of a: vectors y with y.x = 0 for all x in a
    def ann(U):
        return basis([y for y in range(1, 1 << h) if all((y & x).bit_count() % 2 == 0 for x in U)])
    return ann(basis(list(ann(a)) + list(ann(b))))


def dp(lo, hi, s, e):
    """min sum f(increments) over nondecreasing d_i in [lo_i, hi_i], d_0 = s, final = e. Returns (cost, climbs)."""
    INF = (float('inf'), 10 ** 9)
    best = {s: (0.0, 0)}
    for a, b in list(zip(lo, hi)) + [(e, e)]:
        nb = {}
        for d0, (c0, k0) in best.items():
            for d in range(max(a, d0), b + 1):
                c = (c0 + f(d - d0), k0 + (d > d0))
                if c < nb.get(d, INF):
                    nb[d] = c
        best = nb
        if not best:
            return None
    return best[e]


def dp_span(lo, hi):
    """touch-span version: start and end free (cost only between first and last touch)."""
    INF = (float('inf'), 10 ** 9)
    best = {d: (0.0, 0) for d in range(lo[0], hi[0] + 1)}
    for a, b in list(zip(lo, hi))[1:]:
        nb = {}
        for d0, (c0, k0) in best.items():
            for d in range(max(a, d0), b + 1):
                c = (c0 + f(d - d0), k0 + (d > d0))
                if c < nb.get(d, INF):
                    nb[d] = c
        best = nb
    return min(best.values())


def actual(dims, s, e):
    cur, c, k, P = s, 0.0, 0, Counter()
    for d in dims + [e]:
        if d > cur:
            c += f(d - cur); k += 1; P[d - cur] += 1
            cur = d
    return c, k, P


def main():
    Wd = Word(sys.argv[1])
    w = Wd.w
    v, h = Wd.v, Wd.h
    FULL = Wd.full
    restored = {int(q) for q in w.get('helper_final', {})}
    ports = w['ports']
    centre_U = {q: basis(U) for q, U in zip(w['centers'], w['center_frames'])}
    # warp episodes: a register is 'warped' from its first warp on (static pins of it are dropped)
    warped_from = {}
    for t, (sec, e) in enumerate(Wd.events):
        if e[0] == 'warp':
            warped_from.setdefault(e[1], t)
    lend_t = {}
    for t, (sec, e) in enumerate(Wd.events):
        if sec == 'Bwarp':
            lend_t.setdefault(e[1], t)
    cut_t = min(t for t, (sec, e) in enumerate(Wd.events) if sec == 'cut')

    def pins_of(q, t):
        """static (lower, upper) of partner q at time t"""
        if q in warped_from and t >= warped_from[q]:
            return (), FULL
        lo = Wd.start[q]
        if q < v:
            lo = basis(list(lo) + [ports[q]])
        return lo, Wd.end[q]

    rows = []
    for r in Wd.helpers:
        tt = [t for t, s in Wd.touches.get(r, []) if t in Wd.gateframe]
        parts = [('', tt)]
        if r in lend_t:
            parts = [('pre', [t for t in tt if t < lend_t[r]]), ('post', [t for t in tt if t > lend_t[r]])]
        for part, ts in parts:
            if not ts:
                continue
            S = Wd.start[r] if part != 'post' else ()
            E = Wd.end[r] if part != 'pre' else FULL
            if r in restored and part == 'pre':
                E = FULL
            seq = []   # (time, lower, upper, gate frame)
            for t in ts:
                sec, e = Wd.events[t]
                lo, hi = (), FULL
                for q in regs(e):
                    if q == r:
                        continue
                    a, b = pins_of(q, t)
                    lo = basis(list(lo) + list(a))
                    hi = meet(hi, b) if b != FULL else hi
                seq.append((t, lo, hi, Wd.gateframe[t]))
            if r in centre_U and part != 'post':
                U = centre_U[r]
                seq.append((cut_t - 0.5, U, U, U))
                seq.sort(key=lambda x: x[0])
            k = len(seq)
            lo_acc, los = basis(list(S)), []
            for _, a, _, _ in seq:
                lo_acc = basis(list(lo_acc) + list(a))
                los.append(lo_acc)
            hi_acc, his = E, [None] * k
            for i in range(k - 1, -1, -1):
                hi_acc = meet(hi_acc, seq[i][2]) if seq[i][2] != FULL else hi_acc
                his[i] = hi_acc
            feas = all(contained(a, b) for a, b in zip(los, his))
            # positive control: the actual gate frame satisfies the pins
            ctrl = all(contained(lo, F) and contained(F, hi) for (_, _, _, F), lo, hi in zip(seq, los, his))
            dims = []
            for _, _, _, F in seq:
                if not dims or len(F) != dims[-1]:
                    dims.append(len(F))
            ac, ak, aP = actual(dims, len(S), len(E))
            lb = dp([len(x) for x in los], [len(x) for x in his], len(S), len(E))
            lbs = dp_span([len(x) for x in los], [len(x) for x in his])
            span_actual = actual(dims[1:], dims[0], dims[0])[0] if dims else 0.0
            # 8D sub-lemma census
            reads8 = sorted({ports[e[1] - v] for t in ts for (sec, e) in [Wd.events[t]]
                             if e[0] == 'add' and v <= e[1] < 2 * v and e[2] == r and len(Wd.gateframe[t]) == 8})
            reads = sorted({ports[e[1] - v] for t in ts for (sec, e) in [Wd.events[t]]
                            if e[0] == 'add' and v <= e[1] < 2 * v and e[2] == r})
            rows.append(dict(reg=r, part=part, width=len(E) - len(S), actual=ac, actual_climbs=ak,
                             lb=lb[0] if lb else None, lb_climbs=lb[1] if lb else None, feasible=feas,
                             span_actual=span_actual, span_lb=lbs[0], span_lb_climbs=lbs[1],
                             ctrl=ctrl, dims=dims, lo_dims=[len(x) for x in los], hi_dims=[len(x) for x in his],
                             target_ports=len(reads), target_ports_8d=len(reads8), sdim=len(S), edim=len(E)))
    tot_a = sum(x['actual'] for x in rows)
    tot_lb = sum(x['lb'] for x in rows)
    out = dict(word=sys.argv[1], streams=len(rows), all_feasible=all(x['feasible'] for x in rows),
               all_ctrl=all(x['ctrl'] for x in rows),
               actual_logmass=tot_a, lb_logmass=tot_lb, ratio=tot_lb / tot_a,
               actual_climbs=sum(x['actual_climbs'] for x in rows), lb_climbs=sum(x['lb_climbs'] for x in rows),
               streams_lb_equals_actual=sum(1 for x in rows if abs(x['lb'] - x['actual']) < 1e-9),
               streams_lb_zero=sum(1 for x in rows if x['lb'] < 1e-12),
               span_actual=sum(x['span_actual'] for x in rows), span_lb=sum(x['span_lb'] for x in rows),
               span_lb_climbs=sum(x['span_lb_climbs'] for x in rows),
               span_lb_climb_hist=dict(sorted(Counter(x['span_lb_climbs'] for x in rows).items())),
               by_width={})
    byw = defaultdict(lambda: [0, 0.0, 0.0, 0, 0, 0.0, 0.0])
    for x in rows:
        b = byw[x['width']]
        b[0] += 1; b[1] += x['actual']; b[2] += x['lb']; b[3] += x['actual_climbs']; b[4] += x['lb_climbs']
        b[5] += x['span_actual']; b[6] += x['span_lb']
    out['by_width'] = {k: dict(n=b[0], actual=round(b[1], 2), lb=round(b[2], 2), climbs=b[3], lb_climbs=b[4],
                           span_actual=round(b[5], 2), span_lb=round(b[6], 2))
                       for k, b in sorted(byw.items())}
    out['read8d_streams'] = sum(1 for x in rows if x['target_ports_8d'])
    out['read8d_multi_port'] = sum(1 for x in rows if x['target_ports_8d'] > 1)
    out['lb_climb_hist'] = dict(sorted(Counter(x['lb_climbs'] for x in rows).items()))
    out['actual_climb_hist'] = dict(sorted(Counter(x['actual_climbs'] for x in rows).items()))
    print(json.dumps(out, indent=1))
    if '--json' in sys.argv:
        json.dump(rows, open(sys.argv[sys.argv.index('--json') + 1], 'w'))


if __name__ == '__main__':
    main()
