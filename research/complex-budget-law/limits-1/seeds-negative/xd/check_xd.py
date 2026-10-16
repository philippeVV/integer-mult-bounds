"""limits-1/seeds-negative: split Cayley hexagon Tanner codes (expXD), minimal replay.
H(3): 364 + 364, order (3,3), girth 12; natural [4,2] code over F3/F9 has S = 128 (S/E = 0.0879), random gates S <= 3.
H(4) - M (matching seed 0): girth 12; natural code S = 163, one-sided x^2 twist S = 119.
usage: python -I -B check_xd.py   (numpy + numba)"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tanner import *

ok = True
def chk(c, msg):
    global ok; ok &= bool(c); print(('ok   ' if c else 'FAIL ') + msg, flush=True)

H = build(3); NP, NL = incidence(H); g = girth(graph(H), NP + NL)
chk(NP == NL == 364 and g == 12, f'H(3): {NP} points, {NL} lines, girth {g}')
for e in (1, 2):
    T = Tanner(H, GF(3, e), mode='nat'); S = T.dimW()
    chk(T.mds_ok() and S == 128, f'H(3) natural [4,2] over F{3**e}: S = {S}, S/E = {S/T.E:.4f}, a = {a_of(S, T.E):.3f} (b = 11)')
mx = max(Tanner(H, GF(3, 2), mode='randgate', seed=s).dimW() for s in range(3))
chk(mx <= 3, f'H(3) random MDS gates over F9: S <= {mx}')
H4 = build(4); NP4, NL4 = incidence(H4); M = matching(H4, random.Random(0))
g4 = girth(graph(H4, M), NP4 + NL4)
S0 = Tanner(H4, GF(2, 2), drop=M, sp=0, sl=0, mode='nat').dimW()
S1 = Tanner(H4, GF(2, 2), drop=M, sp=0, sl=1, mode='nat').dimW()
chk(g4 == 12 and S0 == 163 and S1 == 119, f'H(4)-M: girth {g4}, natural S = {S0} (S/E = {S0/5460:.4f}), x^2 twist S = {S1} (the twist hurts)')
need = 1 / (1 + 11 ** 0.5)
chk(128 / 1456 < 0.15 and S0 / 5460 < 0.15, f'stop test (S/E < 0.15) fires; a >= 1/2 needs S/E >= {need:.3f}')
print('PASS_EXPXD' if ok else 'FAIL')
