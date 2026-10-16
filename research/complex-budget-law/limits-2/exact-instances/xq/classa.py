"""SA over block (quotient) layouts with per-class freedom: classes i mod n'; state = class permutation, global line offset r,
per-class (m_x, cc_x, lo_x) cyclic sub-layouts.  Objective OBJ = F (true schedule size, min-cost flow) or cost (Boot).
usage: classa.py q n' secs T0 T1 seed OBJ out"""
import sys, math, random, time, numpy as np, importlib.util
q,n1,secs,T0,T1,seed,OBJ,out=int(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),int(sys.argv[6]),sys.argv[7],sys.argv[8]
sys.argv=[sys.argv[0],f'graph{q}.txt']
spec=importlib.util.spec_from_file_location('mcf','mcf.py'); M=importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
N=q*q+q+1; b=N//n1; V=2*N; rng=random.Random(seed); t=np.arange(b)
ms=[m for m in range(1,b) if math.gcd(m,b)==1]
st=dict(perm=list(range(n1)),r=0,m=[1]*n1,cc=[b//3]*n1,lo=[-1]*n1)
def build(st):
    key=np.empty(V,dtype=np.int64)
    for T,x in enumerate(st['perm']):
        key[x+n1*t]=T*10*b+2*((st['m'][x]*t)%b)
        key[N+(x+st['r']+n1*t)%N]=T*10*b+2*((st['m'][x]*t-st['cc'][x])%b)+st['lo'][x]
    return np.argsort(key,kind='stable')
def obj(st):
    c,F=M.F(build(st)); return (F if OBJ=='F' else c),c,F
cur=obj(st); best=cur[0]; print('start',cur,flush=True)
t0=time.time(); it=0
while time.time()-t0<secs:
    it+=1; T=T0*(T1/T0)**((time.time()-t0)/secs)
    new={k_:(v[:] if isinstance(v,list) else v) for k_,v in st.items()}
    mv=rng.random()
    if mv<0.25: i,j=rng.sample(range(n1),2); new['perm'][i],new['perm'][j]=new['perm'][j],new['perm'][i]
    elif mv<0.35: i=rng.randrange(n1); x=new['perm'].pop(i); new['perm'].insert(rng.randrange(n1),x)
    elif mv<0.45: new['r']=rng.randrange(n1)
    else:
        x=rng.randrange(n1); w=rng.random()
        if w<0.3: new['m'][x]=rng.choice(ms)
        elif w<0.8: new['cc'][x]=(new['cc'][x]+rng.choice((-2,-1,1,2)))%b
        else: new['lo'][x]=-new['lo'][x]
    o=obj(new)
    if o[0]<=cur[0] or rng.random()<math.exp(-(o[0]-cur[0])/T):
        st,cur=new,o
        if cur[0]<best:
            best=cur[0]; open(out,'w').write('\n'.join(map(str,build(st)))+'\n'); print(f'{time.time()-t0:.0f}s it {it} best {cur} T {T:.2f}',flush=True)
print('final best',best,'iters',it)
