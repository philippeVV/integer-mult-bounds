"""Price a plane-port decoder with the repo ceiling tool (first-order: non-centre local word unchanged)."""
import sys
from fractions import Fraction as Q
sys.path.insert(0,'./research/complex-local-word-ceiling')
import ceiling as C
WORD={1:2412,2:418,3:432,4:142,5:82,8:8}
def price(dec_children, unpaid=2):
    local={r:n for r,n in WORD.items() if r!=8}
    for r,k in dec_children.items(): local[r]=local.get(r,0)+k
    cen=sum(r*k for r,k in dec_children.items())
    return cen, C.ceiling(local,73,896,9,cen,grid=Q(1,10**12),unpaid=unpaid)[1]
if __name__=='__main__':
    base=C.ceiling(WORD,73,896,9,64,grid=Q(1,10**12))[1]
    print('baseline a_cap',float(base))
    for dc in [{8:8},{9:16},{8:16},{9:16,2:8},{9:18},{2:112},{7:16,2:20}]:
        try:
            c,a=price(dc); print(dc,'C=',c,'a_cap',float(a),'x',float(a/base))
        except ValueError as e: print(dc,'C=',sum(r*k for r,k in dc.items()),'REJECT',e)
