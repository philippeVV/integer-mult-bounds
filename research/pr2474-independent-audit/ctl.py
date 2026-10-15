import os
import json,gzip,sys,copy,io,contextlib,tempfile
from pathlib import Path
PK=Path(os.environ['PKG']+'/new/borrow')
sys.path.insert(0,str(PK));import audit_local as al,word_patch as WP
sys.path.insert(0,'.');import indep
w0=json.loads(gzip.decompress((PK/'words/BORROW-WORD-43base.json.gz').read_bytes()))
tmp=Path(tempfile.mkdtemp())
def pkg(m):
    p=tmp/'m.json.gz';p.write_bytes(gzip.compress(json.dumps(m).encode()))
    try:
        with contextlib.redirect_stdout(io.StringIO()):al.run(p,tmp/'o.json',None,None)
        return 'ACCEPT'
    except Exception as e:return 'reject: '+str(e)[:110]
def mine(m):
    try:indep.check(m);return 'ACCEPT'
    except Exception as e:return 'reject: '+repr(e)[:110]
def rederive(m):
    m=copy.deepcopy(m);A=m['A'];npre=m['prelude_events'];m['A']=A[npre:];m['cleanup']=WP.derive_cleanup(m,WP.restored_of(w0));m['A']=A;return m
ii=1044
cases={}
m=copy.deepcopy(w0);m['B'][ii][3]*=-1;cases['sign flip (pkg-style, cleanup stale)']=m
cases['sign flip, cleanup rederived']=rederive(m)
m=copy.deepcopy(w0);m['B'][ii][3]*=2;cases['doubled (stale)']=m;cases['doubled rederived']=rederive(m)
m=copy.deepcopy(w0);m['warps']=m['warps'][:-1];cases['warp dropped (pkg-style)']=m
# warp dropped but frames otherwise consistent? replace with identity-like: drop warp and rederive cleanup
cases['warp dropped rederived']=rederive(m)
# meaningful semantic mutants
m=copy.deepcopy(w0);m['B'][1630][3]*=-1;cases['sign flip of write 31+=634/4 rederived']=rederive(m)
m=copy.deepcopy(w0);del m['B'][1632]
for r in m['warps']:
    if r[0]>1632:r[0]-=1
cases['last 31 read dropped rederived']=rederive(m)
# move warp earlier than last original use of 31 (B574): pos 500
m=copy.deepcopy(w0);m['warps'][-1][0]=500;cases['warp moved before retirement rederived']=rederive(m)
for k,m in cases.items():print(f'{k:45s} | pkg: {pkg(m)} | mine: {mine(m)}',flush=True)
