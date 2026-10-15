"""Kappa-worth estimate for a favorably closed QQ16 joint-query controller (planning numbers only)."""
import sys
from fractions import Fraction as Q
sys.path.insert(0, str(__import__('pathlib').Path(__import__('os').environ.get('REPO', __import__('pathlib').Path(__file__).resolve().parents[4])) / 'research/complex-local-word-ceiling'))
import ceiling as C
M, FAM, ST = 73, 896, 9
G = Q(1, 10**12)
def cap(p, c=64, unpaid=1):
    try: return C.ceiling(p, M, FAM, ST, c, grid=G, unpaid=unpaid)[1]
    except (ValueError, AssertionError, Exception) as e: return None
k = lambda a: a/(1+a)
SEL = {1:2412,2:418,3:432,4:142,5:82,8:8}
B = {1:2436,2:418,3:432,4:142,5:82}; F = {1:2420,2:418,3:432,4:142,5:82}
MIX = {r: Q(5*F[r]+4*B[r], 9) for r in B}
T74 = {1:2363,2:632,3:326,4:120,5:33,7:4,8:4}
CLAIM = Q('0.001703676387993676')
BIT = Q('0.001733427039589844234749')
BIT1043 = Q('0.001734176288479865766692')
rows = [('selected word, C=64', SEL, 64), ('QQ16 B/F mix, C=64 (PR2079 conv.)', MIX, 64),
        ('QQ16 B/F mix, C=0 (coverage law)', MIX, 0), ('T74 PR2474 local profile, C=64', T74, 64),
        ('T74 profile, centre path free (drop 8s, C=0)', {r:n for r,n in T74.items() if r!=8}, 0)]
res = {}
for name, p, c in rows:
    a = cap(p, c); res[name] = a
    print(f'{name:48s} a_cap={float(a):.9f} kappa_cap={float(k(a)):.9f} vs claim {float(k(a)/CLAIM-1)*100:+.3f}%')
sel = res['selected word, C=64']
for name in rows[1:3]:
    a = res[name[0]]; mult = a/sel
    print(f'multiplier {name[0]}: x{float(mult):.4f}')
print('bit side c (PR589):', float(BIT), f'{float(BIT/CLAIM-1)*100:+.3f}% over claim; (#1043 {float(BIT1043):.12f}, {float(BIT1043/CLAIM-1)*100:+.3f}%)')
