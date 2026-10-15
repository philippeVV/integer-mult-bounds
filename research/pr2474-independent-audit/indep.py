"""Independent re-check of a BORROW-WORD: randomized non-commutative matrix replay (both orientations), own frame chain, paid profile, widths."""
import json,gzip,sys,random
from collections import Counter
PR=(1<<61)-1;DIM=3
def inv(a):return pow(a%PR,PR-2,PR)
def q(n,d):return n%PR*inv(d)%PR
def rmat(rng):
    while True:
        M=[[rng.randrange(PR) for _ in range(DIM)] for _ in range(DIM)]
        I=minv(M)
        if I:return M,I
def minv(M):
    n=DIM;A=[r[:]+[int(i==j) for j in range(n)] for i,r in enumerate(M)]
    for c in range(n):
        p=next((r for r in range(c,n) if A[r][c]),None)
        if p is None:return None
        A[c],A[p]=A[p],A[c];iv=inv(A[c][c]);A[c]=[x*iv%PR for x in A[c]]
        for r in range(n):
            if r!=c and A[r][c]:f=A[r][c];A[r]=[(x-f*y)%PR for x,y in zip(A[r],A[c])]
    return [r[n:] for r in A]
def mv(M,v):return [sum(M[i][j]*v[j] for j in range(DIM))%PR for i in range(DIM)]
def regs(e):return e[1:3] if e[0]=='add' else e[1:2]
def basis(vs):
    piv={}
    for a in vs:
        while a:
            b=a.bit_length()-1
            if b in piv:a^=piv[b]
            else:piv[b]=a;break
    return piv
def span(vs):
    # canonical: reduced row echelon as frozenset of reduced vectors
    piv=basis(vs)
    for k in sorted(piv):
        for j in piv:
            if j>k and piv[j]>>k&1:piv[j]^=piv[k]
    return frozenset(piv.values())
def dim(S):return len(S)
def sub(a,b):return span(list(a)+list(b))==b
def check(w,seed=1):
    v,h,R=w['report']['v'],w['report']['h'],w['report']['R'];n=2*v+R;assert (v,h)==(112,9)
    removed=set(w['removed']);live=[r for r in range(n) if r not in removed]
    restored=set(map(int,w['helper_final']))
    K=len(w['warps']);rng=random.Random(seed)
    mats={}
    for k in range(1,K+1):
        M,I=rmat(rng);mats[k]=M;mats[-k]=I
        M,I=rmat(rng);mats[k+10**6]=M;mats[-(k+10**6)]=I   # dual symbols independent
    npre=w['prelude_events'];Aord=w['A'][npre:]
    byp={}
    for op,r in enumerate(w['warps'],1):byp.setdefault(r[0],[]).append(op)
    fwdB=[]
    for p in range(len(w['B'])+1):
        for op in byp.get(p,[]):fwdB.append(['warp',w['warps'][op-1][1],op])
        if p<len(w['B']):fwdB.append(w['B'][p])
    fwdB+=w['restore']
    # no section other than A-prelude/B/restore may touch removed helpers; check all
    for sec in ('precondition','compensation','A','cut','B','restore'):
        for e in w[sec]:
            rr=[e[1]] if e[0]=='warp' else regs(e)
            assert not any(r in removed or not 0<=r<n for r in rr),(sec,e)
    corr=w['initial_corrections']
    prog=[('corr',t) for t in corr]+[('g',e) for e in w['precondition']+w['compensation']+w['A']+w['cut']+fwdB]
    recorded=[e for e in Aord+fwdB if e[0]=='warp' or not (v<=e[1]<2*v)]
    res={}
    for dual in (False,True):
        init={r:[rng.randrange(PR) for _ in range(DIM)] for r in live}
        val={r:x[:] for r,x in init.items()}
        for c in restored:val[c+n]=[0]*DIM
        DS=10**6
        def warp(r,op):
            if dual:op=(1 if op>0 else -1)*(abs(op)+DS)
            val[r]=mv(mats[op],val[r])
        def add(t,s,c):
            if dual:t,s,c=s,t,(-c)%PR
            val[t]=[(a+c*b)%PR for a,b in zip(val[t],val[s])]
        def scale(t,c):
            if dual:c=inv(c)
            val[t]=[a*c%PR for a in val[t]]
        def gate(e):
            if e[0]=='warp':warp(e[1],e[2])
            elif e[0]=='add':add(e[1],e[2],q(e[3],e[4]))
            else:scale(e[1],q(e[2],e[3]))
        for kind,e in prog:
            if kind=='corr':
                c=__import__('fractions').Fraction(e['coefficient']);c=q(c.numerator,c.denominator);word=e['operator_word']
                if not dual:
                    x=val[e['source']][:]
                    for op in reversed(word):x=mv(mats[op],x)
                    val[e['target']]=[(a+c*b)%PR for a,b in zip(val[e['target']],x)]
                else:
                    x=val[e['target']][:]
                    for op in word:x=mv(mats[-(1 if op>0 else -1)*(abs(op)+DS)],x)
                    val[e['source']]=[(a-c*b)%PR for a,b in zip(val[e['source']],x)]
            else:gate(e)
        if not dual:
            badT=[t for t in range(v) if val[v+t]!=[(a+b)%PR for a,b in zip(init[v+t],init[t])]]
            badC=[c for c in restored if val[c]!=init[c]]
            assert not badT and not badC,(badT[:5],badC[:5])
        # cleanup (derived independently: inverse-reverse of recorded program at FULL, restored -> buffers)
        buf=lambda r:r+n if r in restored else r
        for e in reversed(recorded):
            if e[0]=='warp':warp(buf(e[1]),-e[2])
            elif e[0]=='add':add(buf(e[1]),buf(e[2]),q(-e[3],e[4]))
            else:scale(buf(e[1]),q(e[3],e[2]))
        bad=[]
        for r in live:
            exp=init[r]
            if v<=r<2*v and not dual:exp=[(a+b)%PR for a,b in zip(init[r],init[r-v])]
            if r<v and dual:exp=[(a-b)%PR for a,b in zip(init[r],init[r+v])]
            if val[r]!=exp:bad.append(r)
        assert not bad,('endpoint',dual,bad[:10])
        assert all(val[c+n]==[0]*DIM for c in restored)
        res['dual' if dual else 'forward']='PASS'
    # cleanup equality with shipped cleanup
    derived=[]
    buf=lambda r:r+n if r in restored else r
    for e in reversed(recorded):
        if e[0]=='warp':derived.append(['warp',buf(e[1]),-e[2]])
        elif e[0]=='add':derived.append(['add',buf(e[1]),buf(e[2]),-e[3],e[4]])
        else:derived.append(['scale',buf(e[1]),e[3],e[2]])
    res['cleanup_matches_shipped']=[g[:-1] for g in w['cleanup']]==derived
    # frames: own chain
    full=span([1<<i for i in range(h)])
    cur={r:frozenset() for r in range(n)}
    for i,p in enumerate(w['ports']):cur[i]=span([p])
    for qq,U in w['start_extra'].items():cur[int(qq)]=span(U)
    start=dict(cur);paid=Counter()
    def climb(r,U):
        U=span(U);assert sub(cur[r],U),('nonmonotone',r)
        d=dim(U)-dim(cur[r])
        if d:paid[d]+=1
        cur[r]=U
    for e in w['A'][:npre]:
        if e[0]=='warp':
            _,r,op,fr,to=e;rec=w['warps'][abs(op)-1];assert rec[1]==r
            assert cur[r]==span(fr);G=rec[2]
            gm=lambda x:__import__('functools').reduce(lambda a,i:a^(G[i] if x>>i&1 else 0),range(h),0)
            if op>0:assert span([gm(x) for x in fr])==span(to)
            else:assert span([gm(x) for x in to])==span(fr)
            cur[r]=span(to)
        else:
            for r in regs(e):climb(r,e[-1])
    for e in Aord:
        for r in regs(e):climb(r,e[-1])
    lev={int(t):span(U) for t,U in w['target_levels'].items()}
    assert all(cur[v+t]==lev.get(t,frozenset()) for t in range(v)),'target levels'
    for c,U in zip(w['centers'],w['center_frames']):
        U=span(U);assert cur[c]==U
        served={lev.get(e[1]-v,frozenset()) for e in w['cut'] if e[2]==c}
        for T in sorted(served,key=lambda x:-len(x)):
            assert sub(T,U)
            if T!=U:paid[dim(U)-dim(T)]+=1
            U=T
    warpinfo=[]
    for e in fwdB:
        if e[0]=='warp':
            rec=w['warps'][e[2]-1];r=rec[1];G=rec[2]
            gm=lambda x:__import__('functools').reduce(lambda a,i:a^(G[i] if x>>i&1 else 0),range(h),0)
            assert span([gm(1<<i) for i in range(h)])==full,'G singular'
            assert cur[r]==span(rec[3]),('warp from',e)
            assert span([gm(x) for x in rec[3]])==span(rec[4]),('G(from)',e)
            cur[r]=span(rec[4])
        else:
            for r in regs(e):climb(r,e[-1])
    end={}
    for r in live:
        if r<v or r>=2*v:end[r]=full
        else:end[r]=span([x for x in range(1,1<<h) if (x&w['ports'][r-v]).bit_count()%2==0])
    for qq,U in w['helper_final'].items():end[int(qq)]=span(U)
    for r in live:climb(r,list(end[r]))
    inv_=Counter(dim(end[r])-dim(start[r]) for r in live if r>=2*v)
    res['paid']=dict(sorted(paid.items()));res['inventory']=dict(sorted(inv_.items()));res['helpers']=sum(inv_.values())
    res['paid_matches_report']=res['paid']=={int(k):x for k,x in w['report']['local_profile'].items()}
    res['inv_matches_report']=res['inventory']=={int(k):x for k,x in w['report']['width_inventory'].items()}
    return res
if __name__=='__main__':
    for p in sys.argv[1:]:
        w=json.loads(gzip.decompress(open(p,'rb').read()));print(p.split('/')[-1],check(w),flush=True)
