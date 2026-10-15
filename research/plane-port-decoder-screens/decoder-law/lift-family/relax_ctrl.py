"""Negative controls for relax.py: each drops one true constraint and must find credit >= 10 (so the engine can see cost 71).
(1) infeasible zero pattern al_f = 0 for f in S minus B (violates 2 al alp = 1 there);
(2) pretend g_f.h'_f may vanish (drop the diagonal identity g_f.h'_f = 2/... );
(3) pretend e9.e9' may vanish."""
import random, relax
from relax import maxcredit
from cliques_lib import random_clique, common
rnd = random.Random(8)
while True:
    S = random_clique(8, rnd)
    if len(common(S)) == 8: break
B = S[:7]
v1 = maxcredit(S, B, ZX=S[1:7] + [S[7]], ZY=[], want=True)   # S[7] not in B: infeasible al = 0
print('(1) infeasible al=0 outside B:', v1)
orig = relax.status
relax.status = lambda x, y, ZX, ZY: 'free' if (x == y and x != 'e9') else orig(x, y, ZX, ZY)
v2 = maxcredit(S, S, want=True); print('(2) drop g_f.h_f != 0:', v2)
relax.status = orig
relax.status = lambda x, y, ZX, ZY: 'z' if (x == 'e9' and y == 'e9') else orig(x, y, ZX, ZY)
v3 = maxcredit(S, B, ZX=B, ZY=[], want=True); print('(3) drop e9.e9\' = 1 (B = S minus 1 line, al = 0 on B):', v3)
relax.status = orig
v4 = maxcredit(S, B, ZX=B, ZY=[], want=True); print('    same pattern, true model:', v4)
ok = v1[0] >= 10 and v2[0] >= 10 and v3[0] >= 10 and v4[0] <= 9
print('CONTROLS', 'PASS' if ok else 'FAIL')
