"""Exact obstruction for low leaver sources in ANY circuit of a k8 leaver cell (all relabellings, all orders, any length,
any coefficients), plus the K-block pairing numbers.  Read-only on the repo; GF(2) frames and Fraction ranks only.

Frame rule of the word (cgx_leaver.add/climb, the native frame chains): a gate t += q*s has a frame F containing the
current frames of t and s, and afterwards BOTH t and s have frame F.  Hence after any gate the two operands carry the
same frame, and a register's frame contains the initial frame of every register joined to it by a time-respecting gate
path (its "cone").

Premises, checked below on the given word(s) for all 8 cells (x_j = leaver sources, y_i = targets, paired as bound):
  (P1) init(x_i) is NOT inside hyper(P[y_i])  -> legality: x_i never enters cone(y_i);
  (P2) init(x_k) is NOT inside the block frame of x_j unless x_k, x_j lie in the same block;
  (P3) the cell must deliver C_ij X_j to y_i with C_ij != 0 for all i != j -> x_j enters cone(y_i) for all j != i,
       and nothing outside the cell touches the cell's x's (they have no reads outside the cell).
Lemma (pigeonhole).  Let g_i be the last cell gate touching y_i and v_i its other operand.  Afterwards cone(v_i) =
cone(y_i) = {x_j : j != i}.  v_i is not a target y_m (it would contain x_m, contradicting P1), so v_i is a leaver source.
If v_i = v_k for i != k, then whichever of g_i, g_k is later joins a register whose cone already contains x_i resp. x_k
to that target: impossible by P1.  So i -> v_i is injective from 6 targets into 6 leaver sources, hence bijective: every
leaver source ends with a cone of five sources from all three blocks, so by P2 none ends inside its block frame.
Consequently in every K-block the cell registers are never in I, I is exactly the four retired complement registers, and
dim L = 0 because their protected rows have rank 4 (computed below; their clean forms are rows of the invertible
completion matrix, so this holds for every completion priority, exchange set and carrier placement).

usage: python3 pigeonhole.py HAMMING CODED-CELLS WORD [WORD ...]   (raw word and/or the cold44 input word)
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pairlib import rd, basis, inside, regs, Word


def hyper(p): return basis([x for x in range(1, 512) if (x & p).bit_count() % 2 == 0])


def main():
    ham, cells = sys.argv[1], json.loads(Path(sys.argv[2]).read_text())
    from fractions import Fraction
    assert all((Fraction(c['coefficients'][i][j]) != 0) == (i != j) for c in cells for i in range(6) for j in range(6)), 'P3 coefficients'
    mods = rd(Path(ham) / 'PACKET-PROGRAMS.json.gz'); proto = rd(Path(ham) / 'PROTOTYPE.json')
    blk = {tuple(m['receipt']['block']): b for b, m in enumerate(mods)}
    Fb = [basis([proto['binary_ports'][s] for s in m['receipt']['source_lines']]) for m in mods]
    for wp in sys.argv[3:]:
        w = rd(wp); seq = w['A'] + w['B'] + [w['cleanup'][i][:-1] + [list(f)] for i, f in w.get('early_cleanup', [])]
        print('WORD', wp)
        for c in cells:
            ph = {}
            for k, bd in c['bindings'].items():
                if bd[0] == 'source':
                    b = blk[tuple(bd[1])]; loc = bd[2]
                    ph[k] = ((mods[b]['receipt']['source_lines'][loc] - 8) if loc < 8 else 224 + 16 * b + loc - 8, b)
                else: ph[k] = (112 + bd[1] - 8, None)
            regset = {r for r, _ in ph.values()}; xs = {ph['x%d' % j][0] for j in range(6)}
            cg = [k for k, e in enumerate(seq) if e[0] == 'add' and e[1] in regset and e[2] in regset]
            first, lastg = cg[0], cg[-1]
            # premises on this word (frames at the start of the cell)
            init = {}
            for k, (r, b) in ph.items():
                f = (w['ports'][r],) if r < 112 else basis(w.get('start_extra', {}).get(str(r), []))
                for e in seq[:first]:
                    if r in regs(e): f = basis(e[-1])
                init[k] = f
            for i in range(6):
                assert not inside(init['x%d' % i], hyper(proto['binary_ports'][c['bindings']['y%d' % i][1]])), ('P1', c['cell'], i)
            for j in range(6):
                for k in range(6):
                    same = ph['x%d' % j][1] == ph['x%d' % k][1]
                    assert same or not inside(init['x%d' % k], Fb[ph['x%d' % j][1]]), ('P2', c['cell'], j, k)
            # no gate outside the cell touches a leaver source after the cell starts (P3, second half)
            assert all(not (set(regs(e)) & xs) for k, e in enumerate(seq) if k >= first and k not in set(cg)), ('outside', c['cell'])
            # final partners in the actual word: a bijection targets -> leaver sources
            part = {}
            for i in range(6):
                y = ph['y%d' % i][0]
                k = max(k for k in cg if y in regs(seq[k])); e = seq[k]
                part['y%d' % i] = next(lab for lab, (r, _) in ph.items() if r == (e[2] if e[1] == y else e[1]))
            assert sorted(part.values()) == ['x%d' % j for j in range(6)], (c['cell'], part)
            # final frames of the leaver sources: never inside their block frame
            fin = {}
            for j in range(6):
                r, b = ph['x%d' % j]; f = init['x%d' % j]
                for e in seq[first:]:
                    if r in regs(e): f = basis(e[-1])
                fin['x%d' % j] = (len(f), inside(f, Fb[b]))
                assert not inside(f, Fb[b]), (c['cell'], j)
            print(' ', c['cell'], 'P1 P2 P3 ok; final partners', part, '; leaver-source final frame ranks', {k: v[0] for k, v in fin.items()})
    W = Word(rd(sys.argv[-1]), ham)
    print('PAIRING on', sys.argv[-1])
    for b in range(14):
        d = W.pairing(b)
        print('  block %2d |I| %d rankC %d dimL %d dimK %d rank(LK) %d I=%s' % (b, len(d['I']), d['rankC'], d['dimL'], d['dimK'], d['r'], d['I']))


if __name__ == '__main__':
    main()
