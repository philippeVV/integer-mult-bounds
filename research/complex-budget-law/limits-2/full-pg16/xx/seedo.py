"""seed file from an arbitrary order file (or 'rand:seed'). usage: seedo.py order|rand:s out"""
import sys, random
from common import *
a,out=sys.argv[1],sys.argv[2]
if a.startswith('rand:'):
    o=list(range(V)); random.Random(int(a[5:])).shuffle(o)
else: o=[int(x) for x in open(a).read().split()][:V]
Fv,t,h,b=F(o,True)
with open(out,'w') as f:
    f.write(' '.join(map(str,o))+'\n'); f.write(f'{len(t)}\n')
    for x,y in zip(t,h): f.write(f'{x} {y}\n')
print('boot cost',boot_cost(o),'F',Fv)
