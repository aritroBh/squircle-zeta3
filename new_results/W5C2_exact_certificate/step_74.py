import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
step01 = r'''"""W5C2 step_01 v3: row trimming per sector at p1=16777213, LeanEchelon (fast sparse GE).
Feed order: engine (length-sorted) -> framew5 -> frame234; stop when D_rat & D_sq in span.
Usage: step_01.py pp|pm|mp|mm. Output: work/step01_{sec}.pkl (S, C, rank, fam_counts)."""
import sys, os, time, pickle
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from w5c2_lib import *

secname = sys.argv[1]
sec = {'pp':(1,1),'pm':(1,-1),'mp':(-1,1),'mm':(-1,-1)}[secname]
p = P1
t0 = time.monotonic()
a,b = sec
eng = pickle.load(open(os.path.join(W5C1,f'sector/rows_engine_{a}_{b}.pkl'),'rb'))
f234 = pickle.load(open(os.path.join(W5C1,f'sector/rows_frame234_{a}_{b}.pkl'),'rb'))
fw5 = pickle.load(open(os.path.join(W5C1,f'sector/rows_framew5_{a}_{b}.pkl'),'rb'))
NC = ncols(sec)
feed = [('engine', i, r) for i,r in enumerate(eng)]
feed += [('framew5', i, tuple(sorted(r.items()))) for i,r in enumerate(fw5)]
feed += [('frame234', i, r) for i,r in enumerate(f234)]
log(f'step_01v3 {secname}: feed={len(feed)} (eng {len(eng)}, fw5 {len(fw5)}, f234 {len(f234)}) cols={NC}')
os.makedirs(os.path.join(HERE,'work'), exist_ok=True)
pickle.dump(feed, open(os.path.join(HERE,f'work/feed_{secname}.pkl'),'wb'), protocol=4)
rrows = []
for fam, i, r in feed:
    rrows.append((np.array([c for c,_ in r], np.int32), np.array([int(a_) % p for _,a_ in r], np.int64)))
log(f'step_01v3 {secname}: converted {time.monotonic()-t0:.0f}s')
Dv = {}
for comp in ('rat','sq'):
    D = load_D(sec, comp)
    Dv[comp] = vec_modp(D, p)
state = {'rat': False, 'sq': False}
def chk(E):
    for comp in ('rat','sq'):
        if not state[comp]:
            cols, vals = Dv[comp]
            buf = E.reduce(cols, vals)
            state[comp] = not buf.any()
    log(f'step_01v3 {secname}: stopcheck rank={E.r} in_span={state} el={time.monotonic()-t0:.0f}s')
    return state['rat'] and state['sq']
E = LeanEchelon(NC, p)
_, _, stopped = E.feed(rrows, label=f'{secname}v3', t0=t0, logevery=4000, stopcheck=(2000, chk))
log(f'step_01v3 {secname}: feed done rank={E.r} stopped_at={stopped} in_span={state} elims={E.elims} time={time.monotonic()-t0:.0f}s')
assert state['rat'] and state['sq'], 'D not in span after full feed'
S = list(E.pivrow); C = [int(c) for c in E.pivcol[:E.r]]
from collections import Counter
famc = Counter(feed[i][0] for i in S)
out = {'sec': secname, 'p': p, 'S': S, 'C': C, 'rank': E.r, 'ncols': NC,
       'nrows': len(feed), 'fed': stopped if stopped else len(feed), 'fam_counts': dict(famc)}
log(f'step_01v3 {secname}: |S|={len(S)} families={dict(famc)} |C|={len(C)}')
pickle.dump(out, open(os.path.join(HERE,f'work/step01_{secname}.pkl'),'wb'), protocol=4)
log(f'step_01v3 {secname}: SAVED work/step01_{secname}.pkl')
'''
open('step_01.py','w').write(step01)
import subprocess, time
out = open('step_01.out','wb')
p = subprocess.Popen(['nice','-n','10','/bin/bash','-c','for s in pp pm mp mm; do <project>/.venv/bin/python step_01.py $s || exit 1; done; echo ALL_SECTORS_DONE'],
                     cwd=os.getcwd(), stdout=out, stderr=subprocess.STDOUT, start_new_session=True)
open('step_01.pid','w').write(str(p.pid))
print('launched all sectors, pid', p.pid)
time.sleep(90)
print(open('step_01.out').read())