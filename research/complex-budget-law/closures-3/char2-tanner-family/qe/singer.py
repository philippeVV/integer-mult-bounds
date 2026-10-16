"""Singer difference set for PG(2,16): n=273, D subset Z_n, |D|=17. point i ~ line j iff (j-i) mod n in D."""
def gf_mul(a,b,poly=(1<<12)|(1<<6)|(1<<4)|(1<<1)|1,deg=12):
    r=0
    while b:
        if b&1: r^=a
        b>>=1; a<<=1
        if a>>deg&1: a^=poly
    return r
def singer(q=16):
    n=q*q+q+1; deg=12
    # alpha = x; check primitive
    pw=[1]
    for _ in range(4095): pw.append(gf_mul(pw[-1],2))
    assert pw[4095]==1 and len(set(pw[:4095]))==4095
    def fr(y,k):  # y^(2^k)
        for _ in range(k): y=gf_mul(y,y)
        return y
    def tr(y): return y^fr(y,4)^fr(y,8)
    D=sorted(i for i in range(n) if tr(pw[i])==0)
    assert len(D)==q+1
    diffs={(a-b)%n for a in D for b in D if a!=b}; assert len(diffs)==n-1
    return n,D,pw
if __name__=="__main__":
    n,D,_=singer(); print(n,D)
