import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import numpy as np, time, sys
from circ import *
for mode in ['same','pair']:
    c = Circuit(mode=mode, seed=0)
    rng = np.random.default_rng(1)
    t=time.time(); A = c.transfer(rng); 
    nz=(A!=0); 
    print(mode, 'W',c.W,'time',round(time.time()-t,1),'density',nz.mean().round(4),'diag nonzero',nz.diagonal().mean(), 'rank',rank_mod(A))
    mutual=(nz & nz.T); np.fill_diagonal(mutual,False)
    print('  pairs with mutual reach (2-cycles):', mutual.sum()//2, 'of', c.W*(c.W-1)//2)
    # 0/1 positive gates: reachability
    B = c.transfer(rng, low=1, high=2)
    print('  reach density',(B!=0).mean().round(4))
