"""Two-class theorem on schedule-aligned partitions of the exact q = 9 instance: S = sessions whose source (head) lies in the
first t gates of the causal order sigma, for t = 23, 46, ..., 161; and S = sessions whose sink (tail) does.  Same test as
twoclass9.py (draw 1; other draws only if draw 1 finds a subfield root).  -> twoclass_cut.log"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import twoclass9 as T
S = T.S; c = S.c; W = S.W
Ts = [T.g.transfer(np.random.default_rng(s)) for s in S.GSEEDS]
res = {}
for side, ends in (('head', c.SRC), ('tail', c.SNK)):
    for t in range(23, 182, 23):
        Sset = [k for k in range(W) if c.POS[ends[k]] < t]
        if not 0 < len(Sset) < W: continue
        r = [T.check(Ts[0], Sset, 15)]
        if r[0]: r += [T.check(G, Sset, 15) for G in Ts[1:]]
        res[f'{side}<{t}'] = dict(size=len(Sset), roots=r)
        print(f"{side} in sigma[:{t:3d}]  |S|={len(Sset):3d}  subfield roots d<=15: {r}", flush=True)
json.dump(res, open(os.path.join(HERE, 'twoclass_cut.json'), 'w'), indent=1)
