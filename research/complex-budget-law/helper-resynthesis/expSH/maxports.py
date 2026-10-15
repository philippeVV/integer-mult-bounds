"""Max #ports in a d-dim subspace L (= ports orthogonal to a (9-d)-dim dual Q), d=8..5, exhaustive via bitmasks."""
import json, numpy as np
P=json.load(open('instA.json'))['ports']; dot=lambda a,b:(a&b).bit_count()&1
M=np.zeros((512,2),dtype=np.uint64)
for q in range(1,512):
    m=0
    for i,p in enumerate(P):
        if dot(p,q)==0: m|=1<<i
    M[q,0]=m&((1<<64)-1); M[q,1]=m>>64
pc=lambda a: np.unpackbits(a.view(np.uint8),axis=-1).sum(-1)
print('d=8 max',pc(M[1:]).max())
best2=0
for a in range(1,512):
    x=M[a]&M[a+1:]; best2=max(best2,pc(x).max() if len(x) else 0)
print('d=7 max',best2)
# d=6: triples a<b<c with c not in span(a,b)
best3=0
for a in range(1,512):
    for b in range(a+1,512):
        x=M[a]&M[b]&M[b+1:]
        idx=np.arange(b+1,512); ok=(idx!=(a^b))
        v=pc(x[ok]).max() if ok.any() else 0
        if v>best3: best3=v
print('d=6 max',best3)
