import sys, numpy as np
from common import *
o=[int(x) for x in open(sys.argv[1]).read().split()][:V]; c=cyc_order(11,91)
pos=np.empty(V,int); pos[o]=np.arange(V); pc=np.empty(V,int); pc[c]=np.arange(V)
print('boot',boot_cost(o),'F',F(o))
print('spearman vs cyclic(11,91)', np.corrcoef(pos,pc)[0,1])
# cut profile
j=(pos[NB]<pos[:,None]).sum(1); prof=np.cumsum(17-2*j[o]); print('max cut',prof.max(),'profile samples',prof[::39].tolist())
print('first 40 types', ''.join('P' if v<N else 'L' for v in o[:40]), ' last 40', ''.join('P' if v<N else 'L' for v in o[-40:]))
