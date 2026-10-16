import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
from dfs import dfs_laminar
import itertools
n=5; free={(0,1),(0,2),(1,2),(0,3),(0,4),(3,4)}
B=[[1 if i==j or ((min(i,j),max(i,j)) not in free) else 0 for j in range(n)] for i in range(n)]
E=[(a,b) for a,b in free]+[(b,a) for a,b in free]
# all minimal non-laminar directed subsets (by size, search subsets of size k ascending)
mins=[]
for k in range(2,13):
    for S in itertools.combinations(E,k):
        if any(set(m)<=set(S) for m in mins): continue
        if dfs_laminar(B,list(S))[0] is None: mins.append(S)
    if mins and k>=len(mins[0])+2: break
for m in mins: print(len(m), sorted(m))
print('PASS_MINBOW' if mins and len(mins[0]) == 9 else 'FAIL_MINBOW')
