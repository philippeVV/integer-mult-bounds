"""Step 1 diagnostic: per-block |I|, dim L, dim K, rank(LK) and the registers in I, on a cold44 input word.

usage: python3 diag1.py WORD.json.gz HAMMING_DIR [HELPER-MAP.json]
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pairlib import *

wpath, ham = sys.argv[1], sys.argv[2]
w = rd(wpath); W = Word(w, ham); v = W.v
mods = W.mods


def kind(x, b):
    if x < v:
        ss = [s - 8 for s in mods[b]['receipt']['source_lines']]
        if x in ss: return 'src[l%d]' % ss.index(x)
        ob = next(c for c in range(14) if x in [s - 8 for s in mods[c]['receipt']['source_lines']])
        return 'src(blk%d,l%d)' % (ob, [s - 8 for s in mods[ob]['receipt']['source_lines']].index(x))
    if 224 <= x < 448:
        ob = (x - 224) // 16; loc = 8 + (x - 224) % 16
        return ('fh[l%d]' % loc) if ob == b else 'fh(blk%d,l%d)' % (ob, loc)
    return 'h%d' % x


tot = 0
for b in range(14):
    d = W.pairing(b)
    tot += d['r']
    blk = tuple(mods[b]['receipt']['block'])
    print('BLOCK', b, blk, 'F', d['F'], '|I|', len(d['I']), 'rankC', d['rankC'], 'dimL', d['dimL'], 'dimK', d['dimK'], 'rank(LK)', d['r'])
    print('   I:', [(x, kind(x, b), len(W.cur[x]), 'E' if x in W.enter else '') for x in d['I']])
print('TOTAL rank', tot)
