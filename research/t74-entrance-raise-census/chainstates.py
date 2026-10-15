import os
import sys,json,gzip,tempfile,pickle
from pathlib import Path
D=Path(os.environ.get('PKG','gmx-t74-borrow43-handoff')+'/new/borrow').resolve();sys.path.insert(0,str(D))
import build_word,word_patch,copy
def refit_apply(word0,patch):
    w=copy.deepcopy(word0)
    for sec,i,fr in patch['frames']: w[sec][i][-1]=fr
    w['report']['local_profile']=patch['local_profile']; return w
states={}
for canon,wt,h in [('gauge44','base','eeb6a236'),('extra','extra','cd417b90')]:
    tmp=Path(tempfile.mkdtemp())/'w.json.gz'
    build_word.main(str(D/'inputs'/f'{canon}.canonical.json.gz'),str(D/'inputs/blocks.json'),str(D/'inputs'/f'witness_43_{wt}.json'),str(tmp),h)
    w=json.loads(gzip.decompress(tmp.read_bytes())); st={'built':w}
    rf=D/'inputs'/f'refit_43{wt}.json'
    if rf.exists(): w=refit_apply(w,json.loads(rf.read_text()))
    w=word_patch.apply_kernel(w,json.loads((D/'inputs'/f'kernel_43{wt}.json').read_text()))
    w=word_patch.apply_handoff(w,json.loads((D/'inputs'/f'handoff_43{wt}.json').read_text())); st['handoff']=w
    w=word_patch.apply_handoff(w,json.loads((D/'inputs'/f'entrance_43{wt}.json').read_text())); st['entrance']=w
    w=word_patch.apply_lowfour(w,json.loads((D/'inputs'/f'lowfour_43{wt}.json').read_text())); st['lowfour']=w
    w=word_patch.apply_lowfour(w,json.loads((D/'inputs'/f'alias637_43{wt}.json').read_text())); st['alias637']=w
    states[wt]=st
pickle.dump(states,open('states.pkl','wb'))
print('ok')
