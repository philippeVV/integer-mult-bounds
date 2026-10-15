"""(A) Liveness / interval-graph bound on helper deletion for the #2474 T74 borrow words.

Streams: every data register (source 0..v-1, target v..2v-1) is live from time 0 (it holds an input);
targets stay live to the end (outputs); sources until their last forward touch. Every helper register
is a stream [first touch, last touch]; registers that lend at a B-phase warp are split into a pre-warp
stream and a post-warp (alias) stream. Paired A-prelude warps (warp, reads, inverse warp) are nested
borrows of a live source and are not split.

Two conventions for a helper's first touch:
  literal  : precondition/compensation dirt reads at the zero frame count (time of the read)
  relaxed  : wrapper dirt reads are assumed relocatable (as the alias line does for its deleted helpers);
             the stream starts at its first A/cut/B/restore touch.
Model scope: streams are not uncomputed before the cleanup (the cleanup is the inverse-reverse at FULL,
so intervals are nested LIFO and two streams can share a physical register iff their forward intervals
are disjoint). Then #physical registers >= omega (max clique), data take 224 of them (all live at t=0),
so helper registers H >= omega - 224 and deletable helpers <= 429 - (omega - 224).

usage: liveness.py WORD.json.gz [--json out.json]
"""
import sys, json, math
from collections import Counter
from wordmodel import Word

WRAP = ('precondition', 'compensation')


def streams(W, relaxed):
    v = W.v
    out = []   # (start, end, kind, reg, part)
    lend_t = {}
    for t, (sec, e) in enumerate(W.events):
        if sec == 'Bwarp':
            lend_t.setdefault(e[1], []).append(t)
    for r in W.live:
        tt = W.touches.get(r, [])
        if r < v:
            last = max([t for t, s in tt], default=0)
            if r in lend_t:
                tw = min(lend_t[r])
                pre = [t for t, s in tt if t < tw]
                post = [t for t, s in tt if t >= tw]
                out.append((0, max(pre, default=0), 'source', r, 'pre'))
                out.append((tw, max(post), 'alias', r, 'post'))
            else:
                out.append((0, last, 'source', r, ''))
        elif r < 2 * v:
            out.append((0, W.T_END, 'target', r, ''))
        else:
            ts = [t for t, s in tt if not (relaxed and s in WRAP)]
            if not ts:
                continue
            if r in lend_t:
                tw = min(lend_t[r])
                pre = [t for t in ts if t < tw]
                post = [t for t in ts if t >= tw]
                if pre:
                    out.append((min(pre), max(pre), 'helper', r, 'pre'))
                out.append((tw, max(post), 'alias', r, 'post'))
            else:
                out.append((min(ts), max(ts), 'helper', r, ''))
    return out


def omega(iv):
    ev = []
    for s, e, *_ in iv:
        ev.append((s, 0))      # open before close at equal time: closed intervals
        ev.append((e, 1))
    ev.sort()
    cur = best = 0
    arg = None
    for t, k in ev:
        if k == 0:
            cur += 1
            if cur > best:
                best, arg = cur, t
        else:
            cur -= 1
    return best, arg


def main():
    W = Word(sys.argv[1])
    res = {'word': sys.argv[1], 'helpers': len(W.helpers), 'T': W.T_END}
    sec_start = {}
    for t, (sec, e) in enumerate(W.events):
        sec_start.setdefault(sec, t)
    res['section_start'] = sec_start
    for mode in ('literal', 'relaxed'):
        iv = streams(W, mode == 'relaxed')
        om, at = omega(iv)
        live_at = [x for x in iv if x[0] <= at <= x[1]]
        kinds = Counter(x[2] for x in live_at)
        nh_streams = sum(1 for x in iv if x[2] in ('helper', 'alias'))
        Hmin = om - 2 * W.v
        res[mode] = dict(omega=om, at_time=at, live_kinds=dict(kinds), helper_streams=nh_streams,
                         helper_registers_min=Hmin, deletable_max=len(W.helpers) - Hmin)
        # live-count profile per section boundary (helpers+aliases, retired sources)
        prof = {}
        for sec, t in sec_start.items():
            prof[sec] = dict(Counter(x[2] for x in iv if x[0] <= t <= x[1]))
        res[mode]['profile_at_section_start'] = prof
        # max over time of (live helper streams - retired sources): same number, check
        best = 0
        ts = sorted({x[0] for x in iv} | {x[1] for x in iv})
        for t in ts:
            lh = sum(1 for x in iv if x[2] in ('helper', 'alias') and x[0] <= t <= x[1])
            rs = sum(1 for x in iv if x[2] == 'source' and x[1] < t) + sum(0 for _ in ())
            # a source split at a warp is 'retired' only between pre-end and alias start; the alias
            # stream is counted as helper, so count retired source registers as sources whose
            # (pre or unsplit) stream ended before t
            best = max(best, lh - rs)
        res[mode]['max_live_helpers_minus_retired_sources'] = best
    json.dump(res, sys.stdout, indent=1)
    print()
    if '--json' in sys.argv:
        json.dump(res, open(sys.argv[sys.argv.index('--json') + 1], 'w'), indent=1)


if __name__ == '__main__':
    main()
