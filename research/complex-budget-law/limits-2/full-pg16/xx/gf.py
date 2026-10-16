"""Small GF(p^e) arithmetic with numpy tables, and exact rank / nullspace."""
import numpy as np
IRR = {(2,1):[1,1],(2,2):[1,1,1],(2,3):[1,1,0,1],(2,4):[1,1,0,0,1],(2,5):[1,0,1,0,0,1],
       (2,6):[1,1,0,0,0,0,1],(3,1):[0,1],(3,2):[1,0,1],(5,2):[2,1,1],(7,1):[0,1],(5,1):[0,1]}
# poly coefficients low->high, monic, e.g. F4: 1+x+x^2

class GF:
    def __init__(s,p,e):
        s.p,s.e,s.q=p,e,p**e
        q=s.q
        def dig(a): return [(a//p**i)%p for i in range(e)]
        def und(v): return sum(int(c)*p**i for i,c in enumerate(v))
        irr=IRR[(p,e)]
        add=np.zeros((q,q),dtype=np.int16); mul=np.zeros((q,q),dtype=np.int16)
        for a in range(q):
            da=dig(a)
            for b in range(q):
                db=dig(b)
                add[a,b]=und([(x+y)%p for x,y in zip(da,db)])
                prod=[0]*(2*e)
                for i in range(e):
                    for j in range(e): prod[i+j]=(prod[i+j]+da[i]*db[j])%p
                for k in range(2*e-1,e-1,-1):
                    c=prod[k]
                    if c:
                        prod[k]=0
                        for i in range(e): prod[k-e+i]=(prod[k-e+i]-c*irr[i])%p
                mul[a,b]=und(prod[:e])
        s.add,s.mul=add,mul
        s.neg=np.array([und([(-x)%p for x in dig(a)]) for a in range(q)],dtype=np.int16)
        s.sub=add[:,s.neg]
        inv=np.zeros(q,dtype=np.int16)
        for a in range(1,q): inv[a]=int(np.where(mul[a]==1)[0][0])
        s.inv=inv
    def pw(s,a,n):
        r=1
        for _ in range(n): r=int(s.mul[r,a])
        return r
    def frob(s,a,j=1):  # a^(p^j)
        return s.pw(a,s.p**j) if j else a
    def rref(s,M):
        """returns (R, pivots) reduced row echelon (rows below rank zero)."""
        M=M.copy().astype(np.int16); r=0; piv=[]
        nr,nc=M.shape
        for c in range(nc):
            if r>=nr: break
            nz=np.nonzero(M[r:,c])[0]
            if len(nz)==0: continue
            i=r+nz[0]
            if i!=r: M[[r,i]]=M[[i,r]]
            M[r]=s.mul[s.inv[M[r,c]],M[r]]
            col=M[:,c].copy(); col[r]=0
            rows=np.nonzero(col)[0]
            if len(rows):
                M[rows]=s.sub[M[rows], s.mul[col[rows][:,None], M[r][None,:]]]
            piv.append(c); r+=1
        return M[:r],piv
    def rank(s,M): return len(s.rref(M)[1])
    def nullspace(s,M):
        """basis of {x: M x = 0} as rows"""
        R,piv=s.rref(M); nc=M.shape[1]
        free=[c for c in range(nc) if c not in set(piv)]
        B=np.zeros((len(free),nc),dtype=np.int16)
        for t,f in enumerate(free):
            B[t,f]=1
            for i,pc in enumerate(piv): B[t,pc]=s.neg[R[i,f]]
        return B
