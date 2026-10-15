import os
import sys; sys.path.insert(0,'.')
from frames import *
class Classes:
    def __init__(self,W):
        self.W=W; ev=W.ev; par=list(range(len(ev)))
        def f(a):
            while par[a]!=a: par[a]=par[par[a]]; a=par[a]
            return a
        for r,ks in W.byreg.items():
            for a,b in zip(ks,ks[1:]):
                if ev[a][0]=='gate' and ev[b][0]=='gate' and ev[a][2]==ev[b][2]: par[f(a)]=f(b)
        self.cls={k:f(k) for k in range(len(ev)) if ev[k][0]=='gate'}
        self.members={}
        for k,c in self.cls.items(): self.members.setdefault(c,[]).append(k)
        self.frame={c:ev[c][2] for c in self.members}
        # per register: sequence of items: ('C',cls) or ('W',k)
        self.seq={}
        for r,ks in W.byreg.items():
            s=[]
            for k in ks:
                it=('C',self.cls[k]) if ev[k][0]=='gate' else ('W',k)
                if not s or s[-1]!=it: s.append(it)
            self.seq[r]=s
        self.regs={c:set() for c in self.members}
        for r,s in self.seq.items():
            for i,(t,c) in enumerate(s):
                if t=='C': self.regs[c].add(r)
        self.ka=next(i for i,e in enumerate(ev) if e[0]=='checkA')
        self.ent=set(map(int,W.w['start_extra']))
        self.npre=W.w.get('prelude_events',0)
    def closure(self,seeds,x):
        """minimal forward closure of lifting classes `seeds` by vector x"""
        W=self.W; L=set(); todo=list(seeds); flags=set()
        while todo:
            c=todo.pop()
            if c in L: continue
            if x in span_set(self.frame[c]): continue
            L.add(c)
            for r in self.regs[c]:
                s=self.seq[r]
                for i,(t,cc) in enumerate(s):
                    if t=='C' and cc==c: break
                if i+1<len(s):
                    t2,d=s[i+1]
                    if t2=='W': flags.add('warp')
                    else: todo.append(d)
                else:
                    if not inspan(x,W.end[r]): flags.add('end')
        # A-end checks for targets/centers, prelude
        for c in L:
            for k in self.members[c]:
                ref=W.ev[k][3]
                if ref[0]=='A' and ref[1]<self.npre: flags.add('prelude')
            for r in self.regs[c]:
                if W.v<=r<2*W.v or r in W.w['centers']:
                    # last class of r before checkA
                    s=[it for it in self.seq[r]]
                    lastA=None
                    for k in W.byreg[r]:
                        if k<self.ka and W.ev[k][0]=='gate': lastA=self.cls[k]
                    if lastA==c: flags.add('A-end')
        return L,flags
    def evaluate(self,L,x):
        """profile delta (signed Counter) and raised helpers"""
        W=self.W; dp=Counter(); raised=[]
        touched=set()
        for c in L: touched|=self.regs[c]
        for r in touched:
            s=self.seq[r]
            # frames chain with lift
            def fr(it,lift):
                return self.frame[it[1]]
            old=[];new=[]
            st=W.start[r]; first=s[0]
            raise_=r in self.ent and first[0]=='C' and first[1] in L
            if raise_: raised.append(r)
            o=[st]; n=[basis(list(st)+[x]) if raise_ else st]
            for t,c in s:
                if t=='W': return None
                o.append(self.frame[c]); n.append(basis(list(self.frame[c])+[x]) if c in L else self.frame[c])
            o.append(W.end[r]); n.append(W.end[r])
            for a,b in zip(o,o[1:]):
                d=len(b)-len(a)
                if d: dp[d]-=1
            for a,b in zip(n,n[1:]):
                d=len(b)-len(a)
                assert d>=0
                if d: dp[d]+=1
        return {k:v for k,v in dp.items() if v},raised
def span_set(B):
    s={0}
    for b in B: s|={y^b for y in s}
    return s
def inspan(x,B): return x in span_set(B)
