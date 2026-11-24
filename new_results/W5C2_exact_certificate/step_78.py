import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
step01b = r'''"""W5C2 step_01b: row trimming per sector, LeanEchelon with plateau skip-ahead + fallback.
Feed engine (length-sorted) until rank plateaus (no new pivot in PLATEAU_K rows and >=MINROWS),
then framew5 + frame234; if D not in span, feed the skipped engine tail. stopcheck throughout.
Usage: step_01b.py pp|pm|mp|mm"""
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
NENG = len(eng); NFW = len(fw5)
log(f'step_01b {secname}: feed={len(feed)} (eng {NENG}, fw5 {NFW}, f234 {len(f234)}) cols={NC}')
os.makedirs(os.path.join(HERE,'work'), exist_ok=True)
pickle.dump(feed, open(os.path.join(HERE,f'work/feed_{secname}.pkl'),'wb'), protocol=4)
def conv(r):
    return (np.array([c for c,_ in r], np.int32), np.array([int(a_) % p for _,a_ in r], np.int64))
Dv = {}
for comp in ('rat','sq'):
    D = load_D(sec, comp)
    Dv[comp] = vec_modp(D, p)
state = {'rat': False, 'sq': False}
E = LeanEchelon(NC, p)
def chk():
    for comp in ('rat','sq'):
        if not state[comp]:
            cols, vals = Dv[comp]
            buf = E.reduce(cols, vals)
            state[comp] = not buf.any()
    log(f'step_01b {secname}: stopcheck rank={E.r} in_span={state} el={time.monotonic()-t0:.0f}s')
    return state['rat'] and state['sq']
PLATEAU_K = 8000; MINROWS = 16000
skipped = None
ri = 0; last_rank_row = 0
done = False
while ri < len(feed):
    fam, fi, r = feed[ri]
    cols, vals = conv(r)
    buf = E.reduce(cols, vals)
    if buf.any():
        E.add_pivot(buf, src=ri); last_rank_row = ri
    ri += 1
    if ri % 4000 == 0:
        nnz = [len(x) for x in E.piv_idx[-2000:]] or [0]
        log(f'{secname}b row {ri}/{len(feed)} pivots={E.r} elims={E.elims/1e6:.1f}e6 nnz_avg={sum(nnz)/len(nnz):.0f} el={time.monotonic()-t0:.0f}s')
    if ri % 2000 == 0 and chk():
        done = True; break
    # plateau skip-ahead: only while still in engine region
    if skipped is None and ri >= MINROWS and ri < NENG and (ri - last_rank_row) >= PLATEAU_K:
        log(f'step_01b {secname}: plateau at row {ri} (rank {E.r}), skipping engine[{ri}:{NENG}] for now')
        skipped = (ri, NENG)
        ri = NENG  # jump to framew5
if not done and not chk():
    # fallback: feed skipped engine tail
    if skipped:
        log(f'step_01b {secname}: D not in span; feeding skipped engine rows {skipped}')
        for rj in range(skipped[0], skipped[1]):
            fam, fi, r = feed[rj]
            cols, vals = conv(r)
            buf = E.reduce(cols, vals)
            if buf.any(): E.add_pivot(buf, src=rj)
            if rj % 4000 == 0:
                log(f'{secname}b-tail row {rj} pivots={E.r} el={time.monotonic()-t0:.0f}s')
            if rj % 2000 == 0 and chk():
                done = True; break
    done = done or chk()
log(f'step_01b {secname}: feed done rank={E.r} in_span={state} skipped={skipped} elims={E.elims} time={time.monotonic()-t0:.0f}s')
assert state['rat'] and state['sq'], 'D not in span'
S = list(E.pivrow); C = [int(c) for c in E.pivcol[:E.r]]
from collections import Counter
famc = Counter(feed[i][0] for i in S)
out = {'sec': secname, 'p': p, 'S': S, 'C': C, 'rank': E.r, 'ncols': NC,
       'nrows': len(feed), 'fed': ri, 'skipped': skipped, 'fam_counts': dict(famc)}
log(f'step_01b {secname}: |S|={len(S)} families={dict(famc)} |C|={len(C)}')
pickle.dump(out, open(os.path.join(HERE,f'work/step01_{secname}.pkl'),'wb'), protocol=4)
log(f'step_01b {secname}: SAVED work/step01_{secname}.pkl')
'''
open('step_01b.py','w').write(step01b)
print('step_01b written')
import subprocess
print(subprocess.run(['tail','-3','step_01.out'], capture_output=True, text=True).stdout)