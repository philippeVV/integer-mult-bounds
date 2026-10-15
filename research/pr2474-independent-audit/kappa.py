import json,sys
from decimal import Decimal as D,getcontext
getcontext().prec=90
local={1:2363,2:632,3:326,4:120,5:33,7:4,8:4}
idle={1:21857,3:7308,5:35,48:25752}
calls=666
H={}
for r,c in local.items():H[r]=H.get(r,0)+c*calls
for r,c in idle.items():H[r]=H.get(r,0)+c
pkg=json.load(open(sys.argv[1]))['profile']
print('global profile recomputed',dict(sorted(H.items())),'matches B.json',{int(k):v for k,v in pkg.items()}==H)
m=73;W=66304
print('total mass',sum(r*c for r,c in H.items()),'W*m?',W*(m-1), 'local mass/call',sum(r*c for r,c in local.items()))
def mom(s):return sum(D(c)*(D(r)/m)**(1-s) for r,c in H.items())/W
lo,hi=D(0),D('0.01')
print('mom(0)',mom(D(0)))
for _ in range(280):
    mid=(lo+hi)/2
    if mom(mid)<1:lo=mid
    else:hi=mid
print('root s*',lo)
C=D(sys.argv[2])
print('package kappa',C,'diff',lo-C)
