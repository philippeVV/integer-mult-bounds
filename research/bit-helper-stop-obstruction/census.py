import gzip, json, sys, importlib.util
from collections import Counter, defaultdict
ROOT = sys.argv[1]
W = json.loads(gzip.decompress(open(ROOT+'/research/r439-abc-semantic-bridge/inputs/candidate.json.gz','rb').read()))
spec = importlib.util.spec_from_file_location('xc', ROOT+'/research/e8-bit-p111-tuned/code/xcheck_p.py')
xc = importlib.util.module_from_spec(spec); spec.loader.exec_module(xc)
O = xc.Word(W)
v, R = W['v'], W['R']; n = 2*v+R
print('v',v,'R',R,'n',n,'frames',len(W['frames']), {k:len(W[k]) for k in 'DABC'}, 'ret', len(W['ret']))
print('ADDs', sum(len(g[2]) for k in 'DABC' for g in W[k]) + sum(len(r) for r in W['scatter']))
# paths with touching events
paths = [[f] for f in W['start']]
touch = defaultdict(list)  # (r, frame index in path) -> list of (phase, idx, role)
for ph in ['D','A','scatter','B','C']:
    if ph == 'scatter': continue
    for i,(f,d,src) in enumerate(W[ph]):
        for r in [d]+src:
            if paths[r][-1] != f:
                assert O.nested(paths[r][-1], f) and O.fdim[f] > O.fdim[paths[r][-1]]
                paths[r].append(f)
            touch[(r,len(paths[r])-1)].append((ph,i,'d' if r==d else 's'))
for r,f in enumerate(W['final']):
    if paths[r][-1]!=f: paths[r].append(f)
kind = lambda r: 'src' if r<v else ('tgt' if r<2*v else 'hlp')
prof = Counter()
for p in paths: prof.update(O.fdim[b]-O.fdim[a] for a,b in zip(p,p[1:]))
print('climb profile', dict(sorted(prof.items())))
# rank-1 stops
c = Counter(); lst=[]
for r,p in enumerate(paths):
    dims=[O.fdim[f] for f in p]
    for k in range(1,len(p)-1):
        if dims[k]-dims[k-1]==1 and dims[k+1]-dims[k]==1:
            c[(kind(r),dims[k-1],dims[k],dims[k+1])]+=1
            if dims[k-1]==0: lst.append((r,k,p[k-1],p[k],p[k+1]))
print('consecutive unit climbs (kind,d0,d1,d2):', sorted(c.items()))
print('0->1->2 stops', len(lst))
for r,k,a,b,cc in lst[:400]:
    print(kind(r), r, 'path', [O.fdim[f] for f in paths[r]], 'stop frame', b, W['frames'][b], 'next', W['frames'][cc], 'touch', touch[(r,k)])
json.dump({'paths':paths}, open(sys.argv[2],'w'))
