"""(A') Schedule-independent version of the liveness bound.

Mazurkiewicz-trace argument: two gates that touch a common register in non-commuting roles keep their
relative order under every re-ordering of the gate list. Roles per register: R (read as add source),
A (add-into), S (scale), W (warp on the lender). Commuting pairs: R-R, A-A, S-S; everything else
conflicts. The transitive closure is the dependency order <.

For a gate g, a stream is FORCED LIVE at g if it has a gate <= g and a gate >= g (a source only needs
a gate >= g, since it holds its input from time 0; targets are always live). In every re-ordering all
forced-live streams are simultaneously live when g runs, so

    #physical registers >= max_g F(g),   helper registers >= max_g F(g) - 224.

Wrapper dirt reads (compensation) are dropped from the DAG and streams in the relaxed convention,
kept in the literal one. Frames are ignored (any frame re-solve is allowed): this only weakens the
lower bound, so the deletion cap stays valid.

usage: forced.py WORD.json.gz [relaxed|literal]
"""
import sys, json
from collections import defaultdict
from wordmodel import Word, regs


def maxclique(adj, streams):
    from ortools.sat.python import cp_model
    K = len(adj)
    md = cp_model.CpModel()
    x = [md.NewBoolVar('x%d' % i) for i in range(K)]
    for i in range(K):
        for j in range(i + 1, K):
            if j not in adj[i]:
                md.AddBoolOr([x[i].Not(), x[j].Not()])
    md.Maximize(sum(x))
    sol = cp_model.CpSolver()
    sol.parameters.max_time_in_seconds = 600
    sol.parameters.num_workers = 2
    st = sol.Solve(md)
    assert st == cp_model.OPTIMAL, sol.StatusName(st)
    c = [i for i in range(K) if sol.Value(x[i])]
    for i in c:
        for j in c:
            assert i == j or j in adj[i]
    return c


def build(W, mode, frames):
    v = W.v
    keep = [t for t, (sec, e) in enumerate(W.events) if not (mode == 'relaxed' and sec == 'compensation')]
    idx = {t: i for i, t in enumerate(keep)}
    N = len(keep)
    # role sequence per register
    seq = defaultdict(list)    # r -> [(i, role)]
    for i, t in enumerate(keep):
        sec, e = W.events[t]
        if e[0] == 'warp':
            seq[e[1]].append((i, 'W'))
        elif e[0] == 'add':
            seq[e[1]].append((i, 'A'))
            seq[e[2]].append((i, 'R'))
        else:
            seq[e[1]].append((i, 'S'))
    succ = defaultdict(set)
    for r, s in seq.items():
        blocks = []
        for i, role in s:
            if blocks and blocks[-1][0] == role and role != 'W':
                blocks[-1][1].append(i)
            else:
                blocks.append((role, [i]))
        for (r1, b1), (r2, b2) in zip(blocks, blocks[1:]):
            for x in b1:
                succ[x].update(b2)
    if frames:
        # fixed gate frames: on each register the frame chain is monotone, so a gate at a strictly
        # smaller frame precedes one at a strictly larger frame (profile-preserving re-orderings).
        from wordmodel import contained
        fr = {}
        for i, t in enumerate(keep):
            sec, e = W.events[t]
            if e[0] != 'warp' and sec not in ('precondition', 'compensation', 'cut'):
                fr[i] = W.gateframe[t]
        nfe = 0
        for r, s in seq.items():
            # split the register's history at its warps (frames before/after a warp are not comparable)
            seg, segs = [], []
            for i, role in s:
                if role == 'W':
                    segs.append(seg); seg = []
                elif i in fr:
                    seg.append(i)
            segs.append(seg)
            for sg in segs:
                levels = sorted({fr[i] for i in sg}, key=len)
                for a in levels:
                    for b in levels:
                        if len(a) < len(b) and contained(a, b):
                            pass
                # chain check and edges between consecutive distinct frames
                by = defaultdict(list)
                for i in sg:
                    by[fr[i]].append(i)
                ch = sorted(by, key=len)
                for a, b in zip(ch, ch[1:]):
                    assert contained(a, b), ('frames on a register not a chain', r, a, b)
                    for x in by[a]:
                        for y in by[b]:
                            assert x < y
                            succ[x].add(y); nfe += 1
        print('frame-order edges', nfe, file=sys.stderr)
    pred = defaultdict(set)
    for x, ys in succ.items():
        for y in ys:
            pred[y].add(x)
    desc = [0] * N
    for i in range(N - 1, -1, -1):
        m = 1 << i
        for y in succ[i]:
            m |= desc[y]
        desc[i] = m
    anc = [0] * N
    for i in range(N):
        m = 1 << i
        for x in pred[i]:
            m |= anc[x]
        anc[i] = m
    # streams as gate masks (same split rule as liveness.py)
    lend = {}
    for t, (sec, e) in enumerate(W.events):
        if sec == 'Bwarp' and t in idx:
            lend.setdefault(e[1], idx[t])
    streams = []   # (kind, reg, mask)
    for r in W.live:
        if v <= r < 2 * v:
            continue
        g = [i for i, _ in seq.get(r, [])]
        if not g:
            continue
        kind = 'source' if r < v else 'helper'
        if r in lend:
            tw = lend[r]
            pre = sum(1 << i for i in g if i < tw)
            post = sum(1 << i for i in g if i >= tw)
            if pre:
                    streams.append((kind, r, pre))
            streams.append(('alias', r, post))
        else:
            streams.append((kind, r, sum(1 << i for i in g)))
    return keep, idx, N, desc, anc, streams


def main():
    W = Word(sys.argv[1])
    mode = sys.argv[2] if len(sys.argv) > 2 else 'relaxed'
    v = W.v
    keep, idx, N, desc, anc, streams = build(W, mode, '--frames' in sys.argv)
    # pairwise forced overlap (Helly): a and b overlap in every schedule iff each has a gate that
    # precedes-or-equals a gate of the other (a source only needs: b has a gate <= some gate of a).
    D = []
    for kind, r, m in streams:
        x = 0
        mm = m
        while mm:
            low = mm & -mm
            x |= desc[low.bit_length() - 1]
            mm ^= low
        D.append(x)
    K = len(streams)
    adj = [set() for _ in range(K)]
    for i in range(K):
        for j in range(i + 1, K):
            ki, kj = streams[i][0], streams[j][0]
            ij = bool(D[i] & streams[j][2])   # some gate of i <= some gate of j
            ji = bool(D[j] & streams[i][2])
            if ki == 'source' and kj == 'source':
                ok = True
            elif ki == 'source':
                ok = ji
            elif kj == 'source':
                ok = ij
            else:
                ok = ij and ji
            if ok:
                adj[i].add(j); adj[j].add(i)
    clique = maxclique(adj, streams)
    best = (0, None, None)
    stepF = []
    for gi in range(N):
        a, d = anc[gi], desc[gi]
        F = 2 * v - v   # targets
        cnt = {'source': 0, 'helper': 0, 'alias': 0}
        for kind, r, m in streams:
            if kind == 'source':
                if m & d:
                    cnt['source'] += 1
            elif (m & a) and (m & d):
                cnt[kind] += 1
        F += sum(cnt.values())
        stepF.append(F)
        if F > best[0]:
            best = (F, gi, dict(cnt))
    F, gi, cnt = best
    res_clique = dict(size=len(clique), kinds={k: sum(1 for c in clique if streams[c][0] == k) for k in ('source', 'helper', 'alias')})
    res = dict(frames='--frames' in sys.argv, pairwise_clique=res_clique, pairwise_helper_registers_min=len(clique) + v - 2 * v,
               pairwise_deletable_max=len(W.helpers) - (len(clique) + v - 2 * v),
               edges=sum(len(a) for a in adj) // 2, word=sys.argv[1], mode=mode, gates=N, streams=len(streams), forced_max=F,
               at_gate=keep[gi], at_section=W.events[keep[gi]][0], forced_kinds=cnt,
               helper_registers_min=F - 2 * v, helpers_now=len(W.helpers),
               deletable_max=len(W.helpers) - (F - 2 * v))
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
