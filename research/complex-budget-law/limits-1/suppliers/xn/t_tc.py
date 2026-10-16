import numpy as np, twoclass as tc
f = tc.pmul(tc.pmul([1,1],[1,1,1]), [12345,1])
print(tc.subfield_root_degrees(f, 6))
# charpoly check: companion-free: M diagonal-ish with eigen 1 and a F4 element
rng=np.random.default_rng(0)
M=rng.integers(0,1<<16,(6,6))
cp=tc.charpoly(M); print(len(cp)-1, cp[-1])
# evaluate cp at a point vs det
x=777; v=0
for i,c in enumerate(cp):
    p=1
    for _ in range(i): p=tc.mul(p,x)
    v^=tc.mul(c,p)
print(v==tc.gdet(M^(np.eye(6,dtype=np.int64)*x)))
