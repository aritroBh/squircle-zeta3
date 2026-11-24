import os, subprocess, time
dst = '<project>/research/08_theory/M4i_w5_noL'
with open(dst+'/progress.log','a') as fh:
    fh.write('step_23-24(export+verify): R5_noL.json 69 terms max|num|=2487169 lcm=301916160, 14 Li5 args; symbol/cobracket PASS; fresh 1620d: -2.1052657e-1596.\n')

script = r'''"""M4i step E: cross-check. With c_L = (pL,qL) = L-coefficient of R5_best,
Z := v_noL - v_best - c_L * vL_relation + pL*e_130 + qL*e_131
should be an exact zero combination (R5_noL - R5_best = c_L*(L-rel RHS - L) + exact-zero).
Checks: (i) EXACT full point symbol of Z == 0; (ii) fresh 1620-digit value of Z ~ 1e-1600.
"""
import json, math, pickle, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
from flint import arb, ctx
import mpmath as mp
import w5_symbol as W
from w5_search import evaluate
import nice_search as NS

ZERO, ONE, Q, R, A, B = W.ZERO, W.ONE, W.Q, W.R, W.A, W.B
def log(s): print(s, flush=True)

EXTRA_Z = {'Li5((-2+2*sqrt2))': -2+2*R, 'Li5((17-12*sqrt2))': B**4,
           'Li5((1/2-1/4*sqrt2))': Q(F(1,2))-R/4, 'Li5((1/2+1/4*sqrt2))': Q(F(1,2))+R/4,
           'Li5((1/2-1/2*sqrt2))': Q(F(1,2))-R/2, 'Li5((1/2-3/8*sqrt2))': Q(F(1,2))-3*R/8,
           'Li5((-16+12*sqrt2))': -16+12*R}
NATOM = 146
vb = pickle.loads((HERE/'state_opt2.pkl').read_bytes())['v']
vn = pickle.loads((HERE/'state_noL.pkl').read_bytes())['v']
vL = pickle.loads((HERE/'lrel_state.pkl').read_bytes())['v']
pL, qL = vb[130], vb[131]
log(f'c_L = ({pL}) + ({qL})*sqrt2')
# c_L * vL under (p+q*sqrt2) multiplication
cLvL = [F(0)]*NATOM
for i in range(65):
    a, b = vL[i], vL[65+i]
    cLvL[i] = pL*a + 2*qL*b; cLvL[65+i] = pL*b + qL*a
for k in range(7):
    a, b = vL[132+k], vL[139+k]
    cLvL[132+k] = pL*a + 2*qL*b; cLvL[139+k] = pL*b + qL*a
Z = [vn[s] - vb[s] - cLvL[s] for s in range(NATOM)]
Z[130] += pL; Z[131] += qL
log(f'Z assembled: {sum(1 for x in Z if x)} nonzero slots, max|num|={max(abs(x.numerator) for x in Z if x)}')

# (i) EXACT SYMBOL of Z == 0
sol = NS.solve(list(NS.S0), 'cross')
names = sol['names']
symmap = {nm: s for nm, s, _ in sol['basis']}
recmap = {nm: rec for nm, _, rec in sol['basis']}
extrasQ = [NS.W.Q(F(e['p']), F(e['q'])) for e in
           [e for e in json.loads((HERE/'pool.json').read_text())['pool']
            if NS.W.Q(F(e['p']), F(e['q'])) not in set(NS.S0)]]
for k in range(7):
    z = extrasQ[k]
    symmap[132+k] = W.li_symbol(5, z)
symL = {}
t = {}
for i in range(65):
    a, b = Z[i], Z[65+i]
    if a or b:
        t = W.add(t, W.scale(symmap[names[i]], Q(a, b)))
for k in range(7):
    a, b = Z[132+k], Z[139+k]
    if a or b:
        t = W.add(t, W.scale(symmap[132+k], Q(a, b)))
t = W.clean(t)
assert not t, f'Z symbol NONZERO: {list(t.items())[:3]}'
log('(i) EXACT full point symbol of Z == 0: PASS')

# (ii) FRESH VALUE of Z at 1620 digits
digits = 1620
ctx.prec = math.ceil((digits+140)*math.log2(10)); mp.mp.dps = digits+120
DD = digits+60
def qnum(p, q):
    return mp.mpf(p.numerator)/p.denominator + mp.sqrt(2)*mp.mpf(q.numerator)/q.denominator
total = mp.mpf(0)
for i in range(65):
    a, b = Z[i], Z[65+i]
    if a or b:
        nm = names[i]
        if nm == 'z5':
            v = mp.mpf(arb(5).zeta().mid().str(DD+100, radius=False))
        else:
            v = mp.mpf(evaluate(recmap[nm], DD).mid().str(DD+100, radius=False))
        total += qnum(a, b)*v
for k in range(7):
    a, b = Z[132+k], Z[139+k]
    if a or b:
        v = mp.mpf(evaluate(('li', 5, extrasQ[k]), DD).mid().str(DD+100, radius=False))
        total += qnum(a, b)*v
a, b = Z[130], Z[131]
if a or b:
    Lv = mp.mpf(json.loads((HERE/'Lchi8_5.json').read_text())['value'])
    total += qnum(a, b)*Lv
log(f'(ii) value(Z) at fresh {digits} digits = {mp.nstr(total, 8)}   (expect ~1e-1600: exact zero combination)')
'''
open(dst+'/nol_cross.py','w').write(script)
PY='<project>/.venv/bin/python'
t0=time.monotonic()
r = subprocess.run(['nice','-n','10',PY,'nol_cross.py'], capture_output=True, text=True, cwd=dst, timeout=1750)
print('elapsed %.1fs rc=%d' % (time.monotonic()-t0, r.returncode))
print(r.stdout)
print('STDERR:', r.stderr[-3000:])
