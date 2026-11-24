import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
step02 = r'''"""W5C2 step_02: second-prime (p2=16777199) confirmation per sector.
(a) rank mod p2 of S rows = |S|; (b) D_rat, D_sq in span(S) mod p2;
(c) rank mod p2 of R_{S,C} = |S| (det != 0 mod p2).
Usage: step_02.py pp|pm|mp|mm"""
import sys, os, time, pickle
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from w5c2_lib import *

secname = sys.argv[1]
sec = {'pp':(1,1),'pm':(1,-1),'mp':(-1,1),'mm':(-1,-1)}[secname]
p = P2
t0 = time.monotonic()
d1 = pickle.load(open(os.path.join(HERE,f'work/step01_{secname}.pkl'),'rb'))
S, C, NC = d1['S'], d1['C'], d1['ncols']
rows = pickle.load(open(os.path.join(HERE,f'work/rows_{secname}.pkl'),'rb'))
rS = [rows[i] for i in S]
log(f'step_02 {secname}: |S|={len(S)} |C|={len(C)} p={p}')
# (a)(b): full-column echelon of S at p2
E = Echelon(NC, p)
rrows = []
for fam, i, r in rS:
    cols = np.array([c for c,_ in r], np.int32)
    vals = np.array([int(a) % p for _,a in r], np.int64)
    rrows.append((cols, vals))
E.feed(rrows, label=f'{secname}-S', t0=t0, logevery=4000)
assert E.rank() == len(S), (E.rank(), len(S))
log(f'step_02 {secname}: (a) rank mod p2 = {E.rank()} == |S| OK')
for comp in ('rat','sq'):
    D = load_D(sec, comp)
    cols, vals = vec_modp(D, p)
    rem = E.reduce(cols.tolist(), vals.tolist())
    log(f'step_02 {secname}: (b) D_{comp} remainder nnz={len(rem)} (0 => in span)')
    assert not rem
# (c): C-restricted echelon at p2
Cset = np.array(sorted(C), np.int32)
cpos = {int(c): k for k, c in enumerate(Cset)}
E2 = Echelon(len(C), p)
rrows2 = []
for fam, i, r in rS:
    cc = [cpos[c] for c,_ in r if c in cpos]
    vv = [int(a) % p for c,a in r if c in cpos]
    if cc:
        rrows2.append((np.array(cc,np.int32), np.array(vv,np.int64)))
    else:
        rrows2.append((np.array([],np.int32), np.array([],np.int64)))
E2.feed(rrows2, label=f'{secname}-SC', t0=t0, logevery=4000)
log(f'step_02 {secname}: (c) rank of R_S,C mod p2 = {E2.rank()} (|S|={len(S)})')
assert E2.rank() == len(S)
out = dict(d1)
out['p2'] = p
out['p2_rankS'] = E.rank(); out['p2_rankSC'] = E2.rank()
pickle.dump(out, open(os.path.join(HERE,f'work/step02_{secname}.pkl'),'wb'), protocol=4)
log(f'step_02 {secname}: ALL CHECKS PASSED time={time.monotonic()-t0:.0f}s SAVED')
'''
open('step_02.py','w').write(step02)
print('step_02.py written')
import subprocess
print(subprocess.run(['tail','-3','step_01.out'], capture_output=True, text=True).stdout)