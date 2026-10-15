import os
# Replays verify.py lines 117-131 (after the CRLF-only PLAN hash mismatch at line 116) on the existing work dir, unchanged otherwise.
import gzip,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
P=Path(os.environ['PKG']);work=Path(os.environ['WORK'])
def rd(p):
 b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
env=dict(os.environ,TMPDIR=str(work/'tmp'),PYTHONDONTWRITEBYTECODE='1',PYTHONUTF8='1',PYTHONINTMAXSTRDIGITS='0',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',AUDIT_FAST_SHARDS='2',FINE_ENGINE_CACHE=str(work/'fine-cache'))
def run(name,script,*args):
 with (work/(name+'.log')).open('w') as f:r=subprocess.run([sys.executable,'-B',str(script),*map(str,args)],cwd=work,env=env,stdout=f,stderr=subprocess.STDOUT)
 assert r.returncode==0,(name,(work/(name+'.log')).read_text()[-3000:]);print('PASS',name,flush=True)
assert rd(work/'data-rebound/PROFILE.json')==rd(P/'new/PROFILE.json');print('PROFILE equal')
nat={k:work/('native-new-'+k) for k in ('base','extra')}
layout=work/'parent/layout';runtime=layout/'packing74';runtime.mkdir();data=runtime/'candidate';data.mkdir()
for n in ['replay.py','freeze_spec.py','audit_mixed.py','independent_core.py','geometry_controls.py','audit_fast.py']:shutil.copy2(P/n,runtime/n)
for n in ['code','construction']:shutil.copytree(P/n,runtime/n)
(runtime/'inputs').mkdir();shutil.copy2(P/'inputs/PR589-EXPECTED.json',runtime/'inputs/PR589-EXPECTED.json')
for n in ['PLAN.json','SPEC.json','PROFILE.json','OWNERS.json.gz','DISPATCH.json','CONSTRUCTION.json','GLOBAL-TEMPLATE.json.gz']:shutil.copy2((P/'new' if n in ('PLAN.json','SPEC.json','PROFILE.json') else P)/n,data/n)
shutil.copytree(nat['extra'],runtime/'native-extra')
run('strict-t74',runtime/'replay.py','--data',data,'--work',work/'strict','--local',nat['base']/'local')
run('geometry-controls',runtime/'geometry_controls.py','--source',layout/'joint1148-round45/source/research/complex-k492-owner-repair/vendor/multicold','--strict',work/'strict','--out',work/'GEOMETRY-CONTROLS.json')
r=rd(work/'strict/RESULT.json');assert r['status']=='PASS_STRICT_MIXED_T74_EXACT_FILL' and r['kappa']==rd(P/'CLAIM.json')['kappa'];assert len(rd(work/'GEOMETRY-CONTROLS.json')['rejected'])==6
print(r['status'],r['kappa_decimal'])
