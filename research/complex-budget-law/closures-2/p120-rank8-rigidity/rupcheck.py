"""Minimal independent RUP checker (deletions ignored, which is sound): verifies every lemma of a DRUP proof by unit
propagation and that the formula + lemmas propagate to a conflict (refutation)."""
import gzip, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
if len(sys.argv) < 3: sys.argv[1:] = [os.path.join(HERE, 'nondeg.cnf.gz'), os.path.join(HERE, 'nondeg.drup.gz')]
def opn(f): return gzip.open(f, 'rt') if f.endswith('.gz') else open(f)
def read_cnf(f):
    cls=[]
    for line in opn(f):
        if line[0] in 'cp': continue
        l=list(map(int,line.split()));
        if l: cls.append(l[:-1])
    return cls
cls=read_cnf(sys.argv[1]); lemmas=[]
for line in opn(sys.argv[2]):
    t=line.split()
    if not t or t[0]=='d': continue
    lemmas.append(list(map(int,t[:-1])) if t[-1]=='0' else list(map(int,t)))
nv=max(max(abs(x) for x in c) for c in cls+[l for l in lemmas if l])+1
val=[0]*(2*nv+2)       # val[lit index]: 1 true, -1 false, 0 unassigned
def idx(l): return 2*abs(l)+(l<0)
watch={}
DB=[]
def value(l): v=val[2*abs(l)]; return v if l>0 else -v
def setlit(l,trail): val[2*abs(l)]=1 if l>0 else -1; trail.append(l)
def add(c):
    DB.append(list(c)); k=len(DB)-1
    for l in c[:2]: watch.setdefault(l,[]).append(k)
def propagate(trail, start):
    i=start
    while i<len(trail):
        l=-trail[i]; i+=1          # literal that became false
        ws=watch.get(l,[]); keep=[]
        j=0
        while j<len(ws):
            k=ws[j]; j+=1; c=DB[k]
            if len(c)==1:
                keep.append(k); 
                if value(c[0])==-1: keep.extend(ws[j:]); watch[l]=keep; return False
                continue
            if c[0]==l: c[0],c[1]=c[1],c[0]
            if value(c[0])==1: keep.append(k); continue
            for t in range(2,len(c)):
                if value(c[t])!=-1:
                    c[1],c[t]=c[t],c[1]; watch.setdefault(c[1],[]).append(k); break
            else:
                keep.append(k)
                if value(c[0])==-1: keep.extend(ws[j:]); watch[l]=keep; return False
                if value(c[0])==0: setlit(c[0],trail)
        watch[l]=keep
    return True
units=[]
for c in cls:
    if len(c)==1: units.append(c[0])
    else: add(c)
top=[]
def top_assign(u):
    if value(u)==-1: return False
    if value(u)==0: setlit(u,top)
    return True
ok=all(top_assign(u) for u in units) and propagate(top,0)
if not ok: print('formula refuted at top level'); print('PASS_RUP'); sys.exit()
done=len(top)
for n_,lem in enumerate(lemmas):
    trail=[]; conflict=False
    for l in lem:
        if value(l)==1: conflict=True; break    # lemma satisfied at top level -> trivially implied
        if value(l)==0: setlit(-l,trail)
    if not conflict: conflict = not propagate(trail,0)
    for l in trail: val[2*abs(l)]=0
    if not conflict: print('lemma',n_,'NOT RUP'); sys.exit(1)
    if len(lem)==0: print('empty clause verified'); print('PASS_RUP',n_+1,'lemmas'); sys.exit()
    if len(lem)==1:
        if not top_assign(lem[0]) or not propagate(top,done): print('top-level conflict after lemma',n_); print('PASS_RUP',n_+1,'lemmas'); sys.exit()
        done=len(top)
    else: add(lem)
# no explicit empty clause: check final propagation
print('all lemmas RUP; final top-level propagate conflict:', not propagate(top,done))
