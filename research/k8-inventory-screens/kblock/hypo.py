"""Hypothetical pairing: in the actual k8 leaver HPE word, add K-block cell registers to I as if their final frames
were inside the block frame, and recompute dim L, dim K, rank(LK) exactly.  Also lists the final frame rank of all 24
factory registers of each K-block (which registers could ever join I).
usage: python3 hypo.py HPE_WORD HAMMING CODED-CELLS"""
import sys, json, itertools
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pairlib import *

w = rd(sys.argv[1]); W = Word(w, sys.argv[2]); cells = json.loads(Path(sys.argv[3]).read_text())
mods = W.mods; blk = {tuple(m['receipt']['block']): b for b, m in enumerate(mods)}
def phys(bd):
    b = blk[tuple(bd[1])]; loc = bd[2]
    return (mods[b]['receipt']['source_lines'][loc] - 8) if loc < 8 else 224 + 16 * b + loc - 8
cellregs = {}
for c in cells:
    for k, bd in c['bindings'].items():
        if bd[0] == 'source': cellregs.setdefault(blk[tuple(bd[1])], []).append((phys(bd), c['cell'] + ':' + k))
for b in (0, 1, 2, 11, 12, 13):
    F = W.block_frame(b); base = W.pairing(b)
    ss = [s - 8 for s in mods[b]['receipt']['source_lines']]; fac = ss + list(range(224 + 16 * b, 240 + 16 * b))
    fr = {loc: len(W.cur[x]) for loc, x in enumerate(fac)}
    cr = sorted(cellregs[b])
    print('BLOCK', b, 'I', base['I'], 'final ranks of 24 factory regs (loc:rank)', fr)
    print('   cell registers', cr)
    sop = [x for x, _ in cr if x in (ss[0], ss[3], 224 + 16 * b + 9, 224 + 16 * b + 13)]
    for name, extra in (('+4 notri-low leaver sources', sop), ('+all 8 cell registers', [x for x, _ in cr])):
        d = W.pairing(b, I=base['I'] + extra)
        print('   hypothetical I %-28s |I| %d rankC %d dimL %d dimK %d rank(LK) %d' % (name, len(d['I']), d['rankC'], d['dimL'], d['dimK'], d['r']))
    best = 0
    for x, lab in cr:
        d = W.pairing(b, I=base['I'] + [x]); best = max(best, d['r'])
        print('   hypothetical single +%d (%s): dimL %d dimK %d rank(LK) %d' % (x, lab, d['dimL'], d['dimK'], d['r']))
