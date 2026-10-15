import sys,os,types,importlib.util
D=os.path.join(os.path.dirname(os.path.abspath(__file__)),os.environ.get('P2453','p2453')+'/research/t74-entrance-raise-prices/vendor/pr2387')
spec=importlib.util.spec_from_file_location('ceiling_lp',os.path.join(D,'ceiling_lp.py'));LP=importlib.util.module_from_spec(spec);spec.loader.exec_module(LP)
def root(inv,p,lo=0.0016,hi=0.0018,it=45):
    LP.CASES['t']=dict(inv={w_:inv.get(w_,0) for w_ in LP.WID},pg=dict(p),pe=dict(p))
    for _ in range(it):
        mid=(lo+hi)/2
        if LP.lp_mu_float('t',mid)[0]<1: lo=mid
        else: hi=mid
    return lo
INV0={4:33,5:8,7:62,8:86,9:240}; P0={1:2363,2:632,3:326,4:120,5:33,7:4,8:4}
if __name__=='__main__':
    r0=root(INV0,P0); print('root',r0)
    inv=dict(INV0); inv[8]-=1; inv[7]+=1; p=dict(P0); p[3]-=1; p[2]+=1
    print('own-only k=1', root(inv,p)-r0)
