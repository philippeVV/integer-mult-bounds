"""write graph16.txt and seed (order + arcs) files for cyclic orders. usage: seed.py m c out"""
import sys
from common import *
if __name__=="__main__":
    with open('graph16.txt','w') as f:
        f.write(f'{V} 17\n'); [f.write(' '.join(map(str,nb[v]))+'\n') for v in range(V)]
    m,c,out=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]
    o=cyc_order(m,c); Fv,t,h,b=F(o,True)
    with open(out,'w') as f:
        f.write(' '.join(map(str,o))+'\n'); f.write(f'{len(t)}\n')
        for x,y in zip(t,h): f.write(f'{x} {y}\n')
    print(m,c,'F',Fv,'back',int(b.sum()))
