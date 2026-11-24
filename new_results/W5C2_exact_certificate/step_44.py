import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
step03 = r'''"""W5C2 step_03: exact solve y.R_{S,C} = D_C over Q per sector via Dixon lifting.
Both RHS (D_rat, D_sq) at once. Resume-capable checkpointing.
Usage: step_03.py pp|pm|mp|mm
Output: work/yexact_{sec}.pkl with y_rat, y_sq (Fraction lists over S), heights, steps."""
import sys, os, time, pickle, math
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction as F
from math import lcm
from w5c2_lib import *
import w5c2_dixon as dx

secname = sys.argv[1]
sec = {'pp':(1,1),'pm':(1,-1),'mp':(-1,1),'mm':(-1,-1)}[secname]
t00 = time.monotonic()
d2 = pickle.load(open(os.path.join(HERE,f'work/step02_{secname}.pkl'),'rb'))
S, C, NC = d2['S'], d2['C'], d2['ncols']
rows = pickle.load(open(os.path.join(HERE,f'work/rows_{secname}.pkl'),'rb'))
r = len(S)
log(f'step_03 {secname}: building A=R_S,C^T ({r}x{r}) and B')
Cset = sorted(C); cpos = {c:k for k,c in enumerate(Cset)}
A = np.zeros((r, r), np.int64)
for j, si in enumerate(S):
    for c, a in rows[si][2]:
        k = cpos.get(c)
        if k is not None:
            A[k, j] = a   # A = R^T : A[col_pos, row_j]
    if j % 4000 == 3999: log(f'step_03 {secname}: fill {j+1}/{r}')
log(f'step_03 {secname}: A filled maxentry={int(np.abs(A).max())} nnz={int((A!=0).sum())} ({time.monotonic()-t00:.0f}s)')
# RHS
Ds = {}
Bs = np.zeros((r, 2), np.int64)
Ls = {}
for jj, comp in enumerate(('rat','sq')):
    D = load_D(sec, comp)
    DC = {c: v for c, v in D.items() if c in cpos}
    L = 1
    for v in DC.values(): L = lcm(L, v.denominator)
    Ls[comp] = L
    for c, v in DC.items():
        Bs[cpos[c], jj] = int(v * L)
log(f'step_03 {secname}: B built Ls={Ls} maxB={int(np.abs(Bs).max())}')
# Dixon
ckpt_path = os.path.join(HERE, f'work/dixon_ckpt_{secname}.pkl')
start_step = 0
if os.path.exists(ckpt_path):
    ck = pickle.load(open(ckpt_path,'rb'))
    log(f'step_03 {secname}: resume from checkpoint step={ck["step"]}')
Y, steps, M = dx.dixon_solve(A, Bs, P1, ckpt_path=ckpt_path)
if Y is None:
    log(f'step_03 {secname}: NOT CONVERGED steps={steps}')
    sys.exit(2)
log(f'step_03 {secname}: CONVERGED steps={steps} M_digits={len(str(M))} total={time.monotonic()-t00:.0f}s')
# scale by 1/L and convert
out = {'sec': secname, 'steps': steps, 'Ls': Ls}
for jj, comp in enumerate(('rat','sq')):
    L = Ls[comp]
    y = [F(Y[i][jj][0], Y[i][jj][1]) / L for i in range(r)]
    out['y_'+comp] = y
    nz = sum(1 for a in y if a)
    hn = max((len(str(abs(a.numerator))) for a in y), default=1)
    hd = max((len(str(a.denominator)) for a in y), default=1)
    log(f'step_03 {secname}: y_{comp} nonzero={nz}/{r} max_num_digits={hn} max_den_digits={hd}')
# modular cross-check at p2 (fresh)
ok = dx.check_mod([[ (F(Y[i][j][0],Y[i][j][1])/Ls[('rat','sq')[j]]) for j in range(2)] for i in range(r)], A, Bs, P2)
# note: check_mod checks A@y == B/L mod p2 i.e. A@(Y/L)... construct Bscaled = B (since y=Y/L, A y = B/L; A Y = B)
log(f'step_03 {secname}: cross-check mod p2 (A.Y==B): {ok}')
pickle.dump(out, open(os.path.join(HERE,f'work/yexact_{secname}.pkl'),'wb'), protocol=4)
log(f'step_03 {secname}: SAVED work/yexact_{secname}.pkl')
'''
open('step_03.py','w').write(step03)
print('written')
import subprocess
print(subprocess.run(['tail','-3','step_01.out'], capture_output=True, text=True).stdout)