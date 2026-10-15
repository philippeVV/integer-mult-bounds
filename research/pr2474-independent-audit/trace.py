import os
import json,gzip,sys
P=os.environ['PKG']+'/new/borrow/'
def L(p):return json.loads(gzip.decompress(open(p,'rb').read()))
def regs(e):return e[1:3] if e[0]=='add' else e[1:2]
for wt in ['base','extra']:
  fin=L(P+f'words/BORROW-WORD-43{wt}.json.gz')
  lp=json.load(open(P+f'inputs/lowfour_43{wt}.json'))['fields']
  ap=json.load(open(P+f'inputs/alias637_43{wt}.json'))['fields']
  pri=dict(fin);[pri.__setitem__(k,lp[k]) for k in ap]
  print('=====',wt,'keys',sorted(fin.keys()))
  print('fields identical between alias patch and final:',all(fin[k]==ap[k] for k in ap))
  # B diffs
  gd=[i for i,(a,b) in enumerate(zip(pri['B'],fin['B'])) if a[:-1]!=b[:-1]]
  fd=[i for i,(a,b) in enumerate(zip(pri['B'],fin['B'])) if a[:-1]==b[:-1] and a[-1]!=b[-1]]
  print('B gate diffs',gd,'frame-only diffs',len(fd),fd[:20])
  for i in fd[:20]:print('  frame',i,pri['B'][i][:-1],pri['B'][i][-1],'->',fin['B'][i][-1])
  print('warps prior',len(pri['warps']),'final',len(fin['warps']),'common prefix equal',fin['warps'][:-1]==pri['warps'])
  nw=fin['warps'][-1];print('new warp pos',nw[0],'lender',nw[1],'G',nw[2],'from',nw[3],'to',nw[4])
  print('warp positions sorted?',[w[0] for w in fin['warps']][-6:])
  # forward program with warps at positions
  npre=fin['prelude_events'];A=fin['A']
  def fwd(w):
    out=[('Apre',i,e) for i,e in enumerate(w['A'][:npre])]+[('A',i,e) for i,e in enumerate(w['A'][npre:])]
    byp={}
    for op,r in enumerate(w['warps'],1):byp.setdefault(r[0],[]).append(op)
    for p in range(len(w['B'])+1):
      for op in byp.get(p,[]):out.append(('W',p,['warp',w['warps'][op-1][1],op]))
      if p<len(w['B']):out.append(('B',p,w['B'][p]))
    out+=[('R',i,e) for i,e in enumerate(w['restore'])]
    return out
  for name,w in (('prior',pri),('final',fin)):
    F=fwd(w)
    for r in (31,637):
      ev=[(k,(s,i),e[:-1] if e[0]!='warp' else e) for k,(s,i,e) in enumerate(F) if (e[0]=='warp' and e[1]==r) or (e[0]!='warp' and r in regs(e))]
      print(name,'reg',r,'events',len(ev))
      for x in ev:print('   ',x, '' if x[2][0]=='warp' else ('W' if x[2][1]==r else 'R'))
  # other sections
  for sec in ('precondition','compensation','cut','restore'):
    print(sec,'uses31',[e for e in fin[sec] if 31 in regs(e)][:5],'uses637',[e for e in fin[sec] if 637 in regs(e)][:3])
  print('corrections touching 31/637',[t for t in fin['initial_corrections'] if t['source'] in(31,637) or t['target'] in (31,637)])
  print('ports[31]',fin['ports'][31],'start_extra 31/637',fin['start_extra'].get('31'),pri['start_extra'].get('637'),'helper_final',fin['helper_final'].get('31'),fin['helper_final'].get('637'))
  print('target_levels',fin.get('target_levels'),'centers',fin['centers'])
