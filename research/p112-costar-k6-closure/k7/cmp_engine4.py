"""engineC (C clique scan) vs engine4 (../engine) on random filter-surviving k=7 (set, composition) pairs:
identical output required; every INFEASIBLE certificate is replayed.  usage: cmp_engine4.py SEED NPAIRS
(run in this directory; needs cliques.pkl and cscan.so; env ORBK as in stage1.py)."""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "engine"))   # engine4, engine3
ORBK = os.environ.get("ORBK", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "engine", "orbk.pkl"))
import pickle, random, time, sys
import ext7, engine4, engineC
import verify_cert as V
ext7.init()
reps = pickle.load(open(ORBK,"rb"))[6]
rnd = random.Random(int(sys.argv[1]))
pairs=[]
while len(pairs)<int(sys.argv[2]):
    B=rnd.choice(reps); x=rnd.choice([x for x in range(120) if x not in B])
    for xx,A,s,d in ext7.filt_fast(B):
        if xx==x and A is not None:
            for c in s: pairs.append((A,c))
t4=tc=0; diff=0
for A,c in pairs:
    t=time.time(); r4=engine4.propagate(list(A),list(c),ztime=20.0); t4+=time.time()-t
    t=time.time(); rc=engineC.propagate(list(A),list(c),ztime=20.0); tc+=time.time()-t
    if r4!=rc: diff+=1; print("DIFF",A,c,r4[0],rc[0])
    if rc[0]=="INFEASIBLE": assert V.replay(list(A),list(c),rc[1],rc[2])
print(len(pairs),"diff",diff,"t4 %.4f tc %.4f per pair"%(t4/len(pairs),tc/len(pairs)))
