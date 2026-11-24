import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
step01 = r'''"""W5C2 step_01: row trimming per sector at p1=16777213.
Feed all sector rows (sorted by length) to sparse echelon; record S (pivot rows) and C
(pivot columns). Then check D_rat, D_sq reduce to 0 (in span) and z5/lchi do NOT.
Usage: step_01.py pp|pm|mp|mm
Outputs: work/step01_{sec}.pkl with S, C, rows ref path, stats; appends progress.log."""
import sys, os, time, pickle, json
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from w5c2_lib import *

secname = sys.argv[1]
sec = {'pp':(1,1),'pm':(1,-1),'mp':(-1,1),'mm':(-1,-1)}[secname]
p = P1
t0 = time.monotonic()
rows = load_sector_rows(sec)
NC = ncols(sec)
log(f'step_01 {secname}: {len(rows)} rows, {NC} cols, p={p}; converting')
# save the exact combined ordering for downstream steps
os.makedirs(os.path.join(HERE,'work'), exist_ok=True)
pickle.dump(rows, open(os.path.join(HERE,f'work/rows_{secname}.pkl'),'wb'),
            protocol=4)
# convert to numpy mod-p arrays
rrows = []
for fam, i, r in rows:
    cols = np.array([c for c,_ in r], np.int32)
    vals = np.array([int(a) % p for _,a in r], np.int64)
    rrows.append((cols, vals))
log(f'step_01 {secname}: converted; feeding echelon')
E = Echelon(NC, p)
S, C = E.feed(rrows, label=f'{secname}', t0=t0)
log(f'step_01 {secname}: echelon done rank={E.rank()} time={time.monotonic()-t0:.0f}s')
# membership of D
out = {'sec': secname, 'p': p, 'S': S, 'C': C, 'rank': E.rank(), 'ncols': NC, 'nrows': len(rows)}
for comp in ('rat','sq'):
    D = load_D(sec, comp)
    cols, vals = vec_modp(D, p)
    rem = E.reduce(cols.tolist(), vals.tolist())
    out[f'D_{comp}_rem'] = len(rem)
    log(f'step_01 {secname}: D_{comp} remainder nnz={len(rem)} (0 => in span)')
    assert not rem, f'D_{comp} NOT IN SPAN in {secname}'
# sanity vectors must NOT be in span
a,b = sec
for nm in ('z5','lchi'):
    sv = pickle.load(open(os.path.join(W5C1,f'sector/sanity_{nm}_{a}_{b}.pkl'),'rb'))
    if not sv:
        out[f'sanity_{nm}'] = 'empty'
        continue
    cols, vals = vec_modp(sv, p)
    rem = E.reduce(cols.tolist(), vals.tolist())
    out[f'sanity_{nm}'] = len(rem)
    log(f'step_01 {secname}: sanity {nm} remainder nnz={len(rem)} (expect >0)')
# family composition of S
from collections import Counter
famc = Counter(rows[i][0] for i in S)
out['fam_counts'] = dict(famc)
log(f'step_01 {secname}: |S|={len(S)} families={dict(famc)} |C|={len(C)} total={time.monotonic()-t0:.0f}s')
pickle.dump(out, open(os.path.join(HERE,f'work/step01_{secname}.pkl'),'wb'), protocol=4)
log(f'step_01 {secname}: SAVED work/step01_{secname}.pkl')
'''
open('step_01.py','w').write(step01)
r = __import__('subprocess').run(['nohup','nice','-n','10','/bin/bash','-c',
    'for s in pp pm mp mm; do <project>/.venv/bin/python step_01.py $s; done > step_01.out 2>&1 &'],
    capture_output=True, text=True, cwd='<project>/research/08_theory/W5C2_exact_certificate')
print('launched', r.returncode, r.stderr)
import time; time.sleep(20)
print(open('<project>/research/08_theory/W5C2_exact_certificate/step_01.out').read())