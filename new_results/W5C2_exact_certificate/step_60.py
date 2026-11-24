import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
step01 = r'''"""W5C2 step_01 v2: row trimming per sector at p1=16777213.
Feed order: framew5 (required) -> engine (length-sorted file order) -> frame234.
Stop as soon as BOTH D_rat and D_sq reduce to 0. Save S (pivot source rows), C (pivot cols).
Usage: step_01.py pp|pm|mp|mm"""
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
feed = [('framew5', i, tuple(sorted(r.items()))) for i,r in enumerate(fw5)]
feed += [('engine', i, r) for i,r in enumerate(eng)]
feed += [('frame234', i, r) for i,r in enumerate(f234)]
log(f'step_01v2 {secname}: feed={len(feed)} (fw5 {len(fw5)}, eng {len(eng)}, f234 {len(f234)}) cols={NC}')
os.makedirs(os.path.join(HERE,'work'), exist_ok=True)
pickle.dump(feed, open(os.path.join(HERE,f'work/feed_{secname}.pkl'),'wb'), protocol=4)
rrows = []
for fam, i, r in feed:
    rrows.append((np.array([c for c,_ in r], np.int32), np.array([int(a_) % p for _,a_ in r], np.int64)))
# D vectors for stopcheck
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
    log(f'step_01v2 {secname}: stopcheck rank={E.rank()} in_span={state}')
    return state['rat'] and state['sq']
E = Echelon(NC, p)
stopped = E.feed(rrows, label=f'{secname}v2', t0=t0, logevery=4000, stopcheck=(2000, chk))
log(f'step_01v2 {secname}: feed done rank={E.rank()} stopped_at={stopped} in_span={state} time={time.monotonic()-t0:.0f}s')
assert state['rat'] and state['sq'], 'D not in span after full feed'
S = list(E.pivrow); C = [int(c) for c in E.pivcol[:E.rank()]]
from collections import Counter
famc = Counter(feed[i][0] for i in S)
out = {'sec': secname, 'p': p, 'S': S, 'C': C, 'rank': E.rank(), 'ncols': NC,
       'nrows': len(feed), 'fed': stopped if stopped else len(feed), 'fam_counts': dict(famc)}
log(f'step_01v2 {secname}: |S|={len(S)} families={dict(famc)} |C|={len(C)}')
pickle.dump(out, open(os.path.join(HERE,f'work/step01_{secname}.pkl'),'wb'), protocol=4)
log(f'step_01v2 {secname}: SAVED work/step01_{secname}.pkl')
'''
open('step_01.py','w').write(step01)
# also patch Echelon.feed to support stopcheck (currently only DenseEchelon has it)
lib = open('w5c2_lib.py').read()
old = """    def feed(self, rows, label='', maxsec=1e18, t0=None, start=0):
        for ri in range(start, len(rows)):
            row = rows[ri]
            cols = [c for c,_ in row]; vals = [to_modp(a) for _,a in row]
            rem = self.reduce(cols, vals)
            if rem: self.add_pivot(rem)
            if (ri+1) % 4000 == 0:
                el = time.monotonic()-t0
                nnz = [len(x[0]) for x in self.pivots[-2000:]] or [0]
                log(f'{label} row {ri+1}/{len(rows)} pivots={self.rank()} nnz_avg={sum(nnz)/len(nnz):.0f} max={max(nnz)} ops={self.ops/1e9:.2f}e9 el={el:.0f}s')
                if el > maxsec: return ri+1
        return None"""
new = """    def feed(self, rows, label='', t0=None, logevery=4000, stopcheck=None):
        for ri in range(len(rows)):
            cols, vals = rows[ri]
            rem = self.reduce(cols.tolist(), vals.tolist())
            if rem: self.add_pivot(rem, src=ri)
            if (ri+1) % logevery == 0:
                el = time.monotonic()-t0 if t0 else 0
                nnz = [len(x[0]) for x in self.pivots[-2000:]] or [0]
                log(f'{label} row {ri+1}/{len(rows)} pivots={self.rank()} nnz_avg={sum(nnz)/len(nnz):.0f} max={max(nnz)} ops={self.ops/1e9:.2f}e9 el={el:.0f}s')
            if stopcheck and (ri+1) % stopcheck[0] == 0:
                if stopcheck[1](self):
                    log(f'{label} early stop at row {ri+1} pivots={self.rank()}')
                    return ri+1
        return None"""
assert old in lib
open('w5c2_lib.py','w').write(lib.replace(old, new))
print('step_01 v2 + lib patched')