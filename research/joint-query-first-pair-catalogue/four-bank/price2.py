import sys, json
from fractions import Fraction as Q
sys.path.insert(0,'../../complex-local-word-ceiling'); import ceiling as C
M,FAM,ST=73,896,9
G=Q(1,10**10)
SEL={1:2412,2:418,3:432,4:142,5:82,8:8}
B={1:2436,2:418,3:432,4:142,5:82}
F={1:2420,2:418,3:432,4:142,5:82}
MIX={r:Q(5*F[r]+4*B[r],9) for r in B}
ABIT=Q(1733427039589844234749,10**24)
def add(p,extra):
    q=dict(p)
    for r,n in extra.items(): q[r]=q.get(r,0)+Q(n)
    return q
def cap(p,centre):
    try: return C.ceiling(p,M,FAM,ST,centre,grid=G)[1]
    except ValueError: return Q(0)
SELC=cap(SEL,64)
out={'selected_ceiling':str(SELC),'bit_a':str(ABIT)}
print('selected ceiling',float(SELC),'bit a',float(ABIT))
def bisect(f,lo,hi,target,it=30):
    # f decreasing in x; find x with f(x)=target
    for _ in range(it):
        mid=(lo+hi)/2
        if f(mid)>target: lo=mid
        else: hi=mid
    return lo
for cname,c in [('C64',64),('C0',0)]:
    for pname,p in [('B',B),('MIX',MIX)]:
        base=cap(p,c)
        print(f'\n[{cname} {pname}] base cap {float(base):.9f}  x{float(base/SELC):.4f} over selected')
        if c==64: continue
        # per-replica extra rank-1 children (replace idle): per call n/9
        f1=lambda n: cap(add(p,{1:Q(n)/9}),c)
        f8c=lambda n: cap(add(p,{8:Q(n)/9}),c+Q(8*n)/9)   # centre-like copies: rank-8 child + total
        for nm,f,hi in [('rank-1 children per replica',f1,Q(60000)),('centre-like rank-8 copies per replica',f8c,Q(600))]:
            nb=bisect(f,Q(0),hi,ABIT); ns=bisect(f,Q(0),hi,SELC)
            print(f'  {nm}: keep cap>=bit a up to {float(nb):.1f}; break-even vs selected ceiling at {float(ns):.1f}')
            out[f'{cname}_{pname}_{nm}']=[float(nb),float(ns)]
# concrete 16-coordinate repair scenarios
scen=[]
for c in (64,0):
    for k,label in [(1,'once per replica'),(4,'once per pair (4/replica)'),(9,'every call')]:
        for nm,ex,dc in [('16 rank-1',{1:16},0),('2 rank-8 (ideal packed)',{8:2},0),('one rank-16',{16:1},0),
                         ('16 centre-like copies (rank-8 + total)',{8:16},128)]:
            p=add(MIX,{r:Q(n*k,9) for r,n in ex.items()})
            a=cap(p,c+Q(dc*k,9))
            scen.append((c,label,nm,float(a),float(a/SELC-1)*100))
            print(f'C={c:2d} {label:26s} {nm:40s} cap={float(a):.9f} vs selected {float(a/SELC-1)*100:+.3f}%  {"above bit" if a>ABIT else ""}')
out['scenarios']=scen
json.dump(out,open('price2.json','w'),indent=1)
