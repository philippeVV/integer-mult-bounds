"""expQM / check_qm.py -- exact checks for the cross-level shared-decoder obstruction (README sec.3). Stdlib only.
Run: python3 -I -B check_qm.py [--repo PATH]   (default: the repository containing this file, or IMB_REPO)     (~10 s)

P1  Kronecker overlap lemma (exact, GF(p)): for subspaces A,P <= F^m and B <= F^m',
        dim(A (x) F^m' + P (x) B) = dim A * m' + dim P * dim B - dim(A cap P) * dim B .
    A shared copy serving an outer centre term (A (x) everything inner) and an inner call's copy (P (x) B, inside
    outer child P) saves exactly dim(A cap P) * dim B over paying both; zero unless the host child P meets A.
P2  #1094 complex word (selected; LOCAL-PROGRAMS of the frontier1094 replay), forward and backward:
    (a) every paid child is single-register; the 8 centre children live on the carrier file C, disjoint from the
        704 L registers that carry all 3486 climbs;
    (b) each centre child moves its carrier between ZERO and an 8-dim frame F_k = 415+k (increment = all of F_k);
    (c) content coincidence: the carrier content c_k is shared with exactly one L register (forward: the copy
        source s_k; backward: c_k is gathered fresh at ZERO and only subtracted into the sink).  No climb of s_k
        transforms c_k on a slot of F_k: climbs INTO F_k precede the last write of s_k before the copy (other
        content), climbs after the copy leave F_k upward (increment outside F_k).
P3  #589 bit word (PR589 P111 tuned word): the same for the 9 retained totals (copy rank = dim of the total's
    frame at the snapshot): no climb of a total lies between its last write and the snapshot.
"""
import sys, os, json, gzip, io, tarfile, random, importlib.util
REPO = sys.argv[sys.argv.index('--repo') + 1] if '--repo' in sys.argv else (os.environ.get('IMB_REPO') or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..')))
P = 10007


def rank_mod(rows):
    rows = [r[:] for r in rows if any(r)]
    rk, ncol = 0, (len(rows[0]) if rows else 0)
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % P), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], P - 2, P)
        rows[rk] = [x * inv % P for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % P:
                f = rows[i][c]
                rows[i] = [(x - f * y) % P for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk


def kron(u, v):
    return [a * b % P for a in u for b in v]


def part1(rng):
    n = 0
    for _ in range(300):
        m, mp = rng.randint(3, 7), rng.randint(2, 6)
        def sub(d, amb):
            return [[rng.randrange(P) for _ in range(amb)] for _ in range(d)]
        # force a controlled intersection A cap P of dimension k
        k = rng.randint(0, 3); da, dp = rng.randint(k, m), rng.randint(k, m)
        common = sub(k, m)
        A = common + sub(max(0, da - k), m); Pp = common + sub(max(0, dp - k), m)
        B = sub(rng.randint(0, mp), mp)
        dA, dP, dB = rank_mod(A), rank_mod(Pp), rank_mod(B) if B else 0
        dAP = dA + dP - rank_mod(A + Pp)
        E = [[int(i == j) for j in range(mp)] for i in range(mp)]
        span = [kron(a, e) for a in A for e in E] + [kron(p, b) for p in Pp for b in B]
        lhs = rank_mod(span) if span else 0
        assert lhs == dA * mp + dP * dB - dAP * dB, (m, mp, dA, dP, dB, dAP, lhs)
        n += 1
    print(f'P1 Kronecker overlap lemma: {n} random exact instances OK (saving = dim(A cap P) * dim B)')


def part2():
    tf = tarfile.open(os.path.join(REPO, 'research/community-frontier1094-audit/selected-replay.tar.gz'))
    prog = json.load(gzip.open(io.BytesIO(tf.extractfile(
        'candidate/strict/global/base/providers/gauge44/LOCAL-PROGRAMS.json.gz').read())))
    for side in ('forward', 'backward'):
        f = prog[side]
        climbs = [(i, ev) for i, ev in enumerate(f) if ev[0] == 'climb']
        cen = [(i, ev) for i, ev in enumerate(f) if ev[0] == 'center']
        assert len(climbs) == 3486 and len(cen) == 8 and all(ev[4] == 8 for _, ev in cen)
        Lregs = {ev[2] for _, ev in climbs}
        assert Lregs == set(range(704))
        # carriers are the C file: every scalar touching a carrier names ('C', k)
        for i, ev in cen:
            k = ev[1]
            Fk = 415 + k
            assert {ev[2], ev[3]} == {0, Fk}, ev                 # ZERO <-> F_k, increment = all of F_k (rank 8)
            # the carrier's lifetime block: from its clear to its next clear
            lo = max(j for j in range(i) if f[j][0] == 'scalar' and f[j][1][0] == 'clear' and f[j][1][1] == ['C', k])
            hi = min(j for j in range(i + 1, len(f)) if f[j][0] == 'scalar' and f[j][1][0] == 'clear' and f[j][1][1] == ['C', k])
            block = f[lo:hi + 1]
            if side == 'forward':
                cp = [ev2 for ev2 in block if ev2[0] == 'scalar' and ev2[1][0] == 'copy']
                assert len(cp) == 1 and cp[0][1][1] == ['C', k] and cp[0][2] == Fk
                s = cp[0][1][2][1]; ci = lo + block.index(cp[0])
                # all later uses of the carrier only READ it (adds into L targets): no other register gets c_k
                assert all(ev2[1][0] == 'add' and ev2[1][2] == ['C', k] for ev2 in block[block.index(cp[0]) + 2:-1]
                           if ev2[0] == 'scalar'), k
                writes = [j for j in range(len(f)) if f[j][0] == 'scalar' and isinstance(f[j][1][1], list)
                          and f[j][1][1] == ['L', s]]
                lastw = max(j for j in writes if j < ci)
                nextw = min([j for j in writes if j > ci] + [len(f)])
                into = [j for j, e in climbs if e[2] == s and e[4] == Fk]
                window = [(j, e) for j, e in climbs if e[2] == s and lastw < j < nextw]
                assert into and all(j < lastw for j in into), (k, into, lastw)       # climb into F_k saw other content
                assert all(e[3] == Fk for j, e in window), (k, window)              # in-window climbs leave F_k upward
                print(f'  {side} centre {k}: carrier C{k} <- copy of L{s} at F_k={Fk}; climbs of L{s} into F_k at '
                      f'{into} < last write {lastw}; same-content climbs {[(e[3], e[4], e[5]) for _, e in window]} start AT F_k')
            else:
                gat = [ev2 for ev2 in block if ev2[0] == 'scalar' and ev2[1][0] == 'add' and ev2[1][1] == ['C', k]]
                out = [ev2 for ev2 in block if ev2[0] == 'scalar' and ev2[1][0] == 'add' and ev2[1][2] == ['C', k]]
                assert len(gat) == 22 and all(g[2] == 0 for g in gat) and len(out) == 1 and out[0][2] == Fk
                s = out[0][1][1][1]
                prior = [j for j in range(lo) if f[j][0] == 'scalar' and isinstance(f[j][1][1], list) and f[j][1][1] == ['L', s]]
                assert prior, k                                                       # sink holds other content
                print(f'  {side} centre {k}: carrier C{k} = fresh sum of 22 panel rows at ZERO, lifted to F_k={Fk}, '
                      f'subtracted into L{s} (which carries {len(prior)} earlier writes): c_k is in no L register')
    print('P2 #1094: centre children are single-register (C file) and no climb transforms carrier content on a '
          'slot of F_k -> A cap P = 0 for every host child outside the carrier subtree')


def part3():
    sp = importlib.util.spec_from_file_location('wlib', os.path.join(REPO, 'research/e8-bit-p111-tuned/tune/wlib.py'))
    wlib = importlib.util.module_from_spec(sp); sp.loader.exec_module(wlib)
    w = wlib.load(os.path.join(REPO, 'research/e8-bit-p111-tuned/tune/bitE8mut-p111-tuned.json.gz'))
    F = wlib.Frames(w); dim = F.dim
    cur = list(w['start']); log = {q: [] for q in w['ret']}; step = 0
    for sec in ('D', 'A'):
        for fr, t, ss in w[sec]:
            step += 1
            for r in [t, *ss]:
                if r in log and cur[r] != fr:
                    log[r].append(('climb', step, cur[r], fr))
                cur[r] = fr
            if t in log:
                log[t].append(('write', step, fr))
    snapdims = [dim[cur[q]] for q in w['ret']]
    assert snapdims == [7] * 9
    for q in w['ret']:
        ev = log[q]; lw = max(i for i, e in enumerate(ev) if e[0] == 'write')
        assert not [e for e in ev[lw + 1:] if e[0] == 'climb'], q
    # after the snapshot every climb of a total starts at or above its snapshot frame (ascent out of F_q)
    snap = {q: cur[q] for q in w['ret']}
    for sec in ('B', 'C'):
        for fr, t, ss in w[sec]:
            for r in [t, *ss]:
                if r in snap:
                    assert F.le(snap[r], fr)
    print('P3 #589: 9 retained totals copied at rank 7; no climb between last write and snapshot; later climbs '
          'ascend out of the snapshot frame -> no same-content slot overlap')


def main():
    part1(random.Random(7))
    part2()
    part3()
    print('PASS_EXPQM_CHECK')


if __name__ == '__main__':
    main()
