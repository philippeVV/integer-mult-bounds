import os
import json,gzip,sys,importlib.util,collections
P=os.environ['PKG']+'/'
def rd(p):
    b=open(P+p,'rb').read();return json.loads(gzip.decompress(b) if p.endswith('.gz') else b)
def span_dim(vs):
    piv={}
    for a in vs:
        while a:
            b=a.bit_length()-1
            if b in piv:a^=piv[b]
            else:piv[b]=a;break
    return len(piv)
# widths from the shipped words, my own derivation: residual width = dim(final) - dim(start_extra) over remaining helpers in register order
W={}
for key,f in (('gauge44','base'),('extra','extra')):
    w=rd(f'new/borrow/words/BORROW-WORD-43{f}.json.gz');rem=[r for r in range(224,704) if r not in set(w['removed'])]
    se={int(k):v for k,v in w['start_extra'].items()};hf={int(k):v for k,v in w['helper_final'].items()}
    W[key]=[(span_dim(hf[r]) if r in hf else 9)-span_dim(se.get(r,[])) for r in rem]
prov=rd('new/BORROW-PROVIDERS.json')
print('my widths == BORROW-PROVIDERS',{k:W[k]==prov[k]['widths'] for k in W},'R',len(W['gauge44']),collections.Counter(W['gauge44']))
plan=rd('inputs/PLAN.json');parent=rd('inputs/GLOBAL-TEMPLATE.json.gz');T=plan['T']
order=rd('new/borrow-fill/ORDER.json');disp={(s,r):p for s,r,p in rd('new/borrow-fill/DISPATCH.json')}
print('order len',len(order),'perm of pinned calls',sorted(map(tuple,order))==sorted((s,r+T*g) for s,r in plan['order'] for g in range(2)),'dispatch',collections.Counter(disp.values()))
sp=importlib.util.spec_from_file_location('bp',P+'vendor/build_parent.py');mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod)
windows,_=mod.actual_windows(parent,order)
fill=rd('new/borrow-fill/FILL.json.gz');rows=fill['rows'];assert not fill['dispatch_overrides']
R=len(W['gauge44'])
percall=collections.defaultdict(list);bins=collections.defaultdict(list)
for t,hlp,rep,ch,port,off in rows:
    call=tuple(order[t]);wd=W[disp[call]][hlp]
    percall[call].append((hlp,rep,ch,port))
    cap,allowed=windows[f'{rep},{ch}']
    assert t in allowed,('window',t,rep,ch)
    assert 0<=port<112 and 0<=off and off+wd<=cap,('cap',)
    bins[rep,ch,port].append((off,off+wd,t,hlp))
bad=0
for call,lst in percall.items():
    hs=[x[0] for x in lst]
    if sorted(hs)!=list(range(R)):bad+=1
    fam=[rep*896+(2*plan['pairs'][str(ch)]+plan['sides'][str(ch)])*112+port for _,rep,ch,port in lst]
    if len(set(fam))!=len(fam):bad+=1
    if len({(rep,ch,port) for _,rep,ch,port in lst})!=len(lst):bad+=1
print('calls',len(percall),'bad calls',bad)
idle=collections.Counter();ov=0
for key,(cap,times) in windows.items():
    rep,ch=map(int,key.split(','))
    for port in range(112):
        iv=sorted(bins.get((rep,ch,port),[]))
        e=0
        for a,b,_,_ in iv:
            if a!=e:ov+=1
            e=b
        if cap-e:idle[cap-e]+=1
print('overlap/gap count',ov,'bins with roles outside windows',len(set(bins)-{(int(k.split(',')[0]),int(k.split(',')[1]),p) for k in windows for p in range(112)}))
print('idle',dict(sorted(idle.items())))
H=collections.Counter(idle)
for call in map(tuple,order):
    for r,n in prov[disp[call]]['actual_paid_profile'].items():H[int(r)]+=n
print('global',dict(sorted(H.items())),'mass',sum(r*n for r,n in H.items()),'live',plan['live'],'deficit check',73*plan['live']*2-sum(r*n for r,n in H.items()),2*plan['deficit'])
print('max idle rank',max(idle),'window caps',collections.Counter(c for c,_ in windows.values()))
