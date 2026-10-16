import sys, numpy as np
sys.argv = ['e8pad.py', sys.argv[1], sys.argv[2], '0']
exec(open('e8pad.py').read().split("rng = random.Random(0)")[0])
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from lib import clique_ls
zmask = (Z != 0).any(1)
for nm, sub in (('z-only', np.nonzero(zmask)[0]), ('all, cold start', np.arange(len(P)))):
    A = adj[np.ix_(sub, sub)]
    if nm == "z-only" and len(sys.argv) > 3: continue
    best = clique_ls(A, iters=30000, restarts=4, seed=1)
    print(nm, 'cands', len(sub), 'best clique', len(best), 'z!=0 in it', int(zmask[sub[best]].sum()), flush=True)
