"""Check the hand-built zero-helper-cost schedule for the F_2^3 toy interface (README section 3, Toy)."""
def span(vs):
    S={0}
    for v in vs: S|={x^v for x in S}
    return frozenset(S)
perp=lambda p: frozenset(x for x in range(8) if (x&p).bit_count()%2==0)
fr={'s1':span([1]),'s2':span([2]),'s4':span([4]),'t1':span([]),'t2':span([]),'t4':span([]),'C':span([]),'H':span([])}
cap={'t1':perp(1),'t2':perp(2),'t4':perp(4),'C':span([1,2])}
reach={k:set() for k in fr}; reach.update(s1={1},s2={2},s4={4})
Hframes=set()
def gate(dst,src,post):
    G=span(list(fr[dst])+list(fr[src]))
    for r in (dst,src):
        assert G<=cap.get(r,frozenset(range(8))),(dst,src)
        if not post: assert not r.startswith('t'), 'target touched before cut'
        fr[r]=G
    reach[dst]|=reach[src]
    if 'H' in (dst,src): Hframes.add(G)
for d,s in [('H','s2'),('C','s1'),('C','s2')]: gate(d,s,False)
assert reach['C']>={1,2} and fr['C']==cap['C']
for d,s in [('t4','s1'),('t4','s2'),('t2','s4'),('t2','s1'),('t1','H'),('t1','s4')]: gate(d,s,True)
need={'t1':{2,4},'t2':{1,4},'t4':{1,2}}
assert all(reach[t]>=need[t] for t in need)
assert len(Hframes)==1
print('PASS toy: interface met, helper H frame-constant at',sorted(next(iter(Hframes))),'-> zero helper cost')
