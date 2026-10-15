"""Price the information-minimal repair (>= 62 coordinates per owner per nine-call replica) against the
selected complex primitive with the pinned ceiling.py.  Per-call profiles are replica averages (fractional
counts allowed, as in the retained model).  Conventions:
  C64: PR2079's convention (centre charge 64 kept although the QQ16 cores are centre-free);
  C0 : coverage law with C = copied centre mass (centre-cut-pair-sharing README), C = 0 for QQ16.
Repair prices per coordinate per live call:
  line : one rank-1 child replacing idle mass (in-place write, no new live data);
  centre: one rank-8 child + 8 units of total mass (same price as a centre coordinate of the selected word);
  family: 72 units of total mass paid by ideal rank-51 idle children (a new live family that is not work)."""
import sys, json
from fractions import Fraction as Q
sys.path.insert(0,'../../complex-local-word-ceiling'); import ceiling as C
M,FAM,ST=73,896,9
G=Q(1,10**10)
SEL={1:2412,2:418,3:432,4:142,5:82,8:8}
B={1:2436,2:418,3:432,4:142,5:82}; F={1:2420,2:418,3:432,4:142,5:82}
MIX={r:Q(5*F[r]+4*B[r],9) for r in B}
ABIT=Q(1733427039589844234749,10**24)
def cap(p,c):
    try: return C.ceiling(p,M,FAM,ST,c,grid=G)[1]
    except ValueError as e: return None
SELC=cap(SEL,64)
def priced(base_c,coords,L,kind):
    n=Q(coords*L,9)  # coordinate-calls per call
    p=dict(MIX); c=Q(base_c)
    if kind=='line': p[1]=p[1]+n
    elif kind=='centre': p[8]=p.get(8,0)+n; c+=8*n
    elif kind=='family': c+=72*n
    return cap(p,c)
rows=[]
for base_c in (64,0):
    for coords in (62,64):
        for kind,Ls in (('line',(1,)),('centre',(1,2,9)),('family',(1,))):
            for L in Ls:
                a=priced(base_c,coords,L,kind)
                s='infeasible/no saving' if a is None else f'{float(a):.9f}  {float(a/SELC-1)*100:+8.3f}% vs selected ceiling  {"> bit a" if a>ABIT else "< bit a"}'
                print(f'C={base_c:2d} coords={coords} {kind:6s} L={L}: {s}')
                rows.append([base_c,coords,kind,L,None if a is None else str(a)])
# break-even per-coordinate-call total-mass price (paid as rank-8 children) for 62 coords, C0
def be(target):
    lo,hi=Q(0),Q(20)
    for _ in range(30):
        mid=(lo+hi)/2
        n=Q(62,9)*mid/8
        p=dict(MIX); p[8]=n
        a=cap(p,8*n)
        if a is not None and a>target: lo=mid
        else: hi=mid
    return lo
b1,b2=be(ABIT),be(SELC)
print('C0, 62 coordinates/replica: max centre-like mass per coordinate-call keeping cap>bit a:',float(b1),' break-even vs selected:',float(b2))
json.dump(dict(selected_ceiling=str(SELC),bit=str(ABIT),rows=rows,breakeven_bit=str(b1),breakeven_sel=str(b2)),open('price3.json','w'),indent=1)
