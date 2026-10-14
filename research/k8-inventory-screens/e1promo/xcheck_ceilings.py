"""Cross-check e1lib backward ceilings at every entering helper's first gate against PR1148's frame_model hi[] (independent
Zassenhaus-meet implementation, downloaded read-only to pr1148/ee/).  (expB)"""
import sys, os, importlib.util
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import e1lib as L
s = importlib.util.spec_from_file_location('frame_model', os.environ.get('FRAME_MODEL', os.path.join(HERE, 'pr1148/ee/frame_model.py')))
M = importlib.util.module_from_spec(s); s.loader.exec_module(M)
w, sha = L.load(sys.argv[1]); W = L.Word(w); W.sha = sha; W.prepare(); W.compute_cones(sys.argv[2]); C, gc = W.ceilings()
m = M.Model(w, relax_given=True)
agree = Counter()
for q in W.entering:
    g = W.revs[q][0]
    hi = m.hi[('G', g)]
    mine = gc[g]
    agree[(L.dim(W.start[q]), L.dim(C[q]), L.mask(hi) == mine, L.mask(hi) == C[q] or L.dim(C[q]) < L.dim(mine))] += 1
print('(entrance dim, my start ceiling dim, first-gate hi equal, consistent):', dict(agree))
# all gates: my gate ceilings vs frame_model hi
eq = sum(L.mask(m.hi[('G', i)]) == gc[i] for i in range(len(W.ev)))
print('gate ceilings equal on', eq, 'of', len(W.ev), 'gates')
