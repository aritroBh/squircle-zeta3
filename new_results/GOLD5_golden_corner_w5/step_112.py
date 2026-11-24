step05 = r'''"""step_05.py -- GOLD5 number level, part 2: LLL fits of the symbol-invisible remainder
R = W5 - N5 - P5 at 1150 and 1400 digits (flint LLL, gram='exact'), over
  17 ring atoms + 10 golden invisible constants + {L(chi5,5), sqrt5*L(chi5,5)}
with sqrt5-multiples of everything (coefficients in Q(sqrt5)).
Acceptance: identical relation vector at both precisions AND residual at working precision.
Controls: planted duplication relation 16Li5(rho)+16Li5(-rho)-Li5(rho^2)=0 must be recovered.
Also: (i) necessity test -- same fit with L-columns dropped; (ii) ZQ5 ladder trade --
replace {L, sqrt5 L} by {Li5(rho^3), Li5(rho^6), their sqrt5 multiples} and refit.
"""
import sys, json, time, mpmath as mp
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
DPS = 1490
mp.mp.dps = DPS + 40
t0 = time.time()
from flint import fmpz_mat, fmpz

def log(msg):
    print(msg, flush=True)
    with open(GOLD + "progress.log", "a") as f:
        f.write(msg + "\n")

nv = {k: mp.mpf(v) for k, v in json.load(open(GOLD + 'number_values.json')).items()}
r5 = mp.sqrt(5)
l2, lf, pi2, pi4 = nv['l2'], nv['lf'], nv['pi']**2, nv['pi']**4
z2, z3, z5 = nv['z2'], nv['z3'], nv['z5']
B = []
for a in range(6):
    B.append((f'l2^{5-a}*lf^{a}', l2**(5-a)*lf**a))
for a in range(4):
    B.append((f'pi2*l2^{3-a}*lf^{a}', pi2*l2**(3-a)*lf**a))
B += [('pi4*l2', pi4*l2), ('pi4*lf', pi4*lf)]
for a in range(3):
    B.append((f'z3*l2^{2-a}*lf^{a}', z3*l2**(2-a)*lf**a))
B += [('z2*z3', z2*z3), ('z5', z5)]
B += [('Li5(rho)', nv['Li5(rho)']), ('Li5(-rho)', nv['Li5(-rho)']), ('Li5(rho^2)', nv['Li5(rho^2)'])]
B += [('K4m*l2', nv['K4m']*l2), ('K4m*lf', nv['K4m']*lf), ('K4q*l2', nv['K4q']*l2), ('K4q*lf', nv['K4q']*lf)]
B += [('K3m*l2^2', nv['K3m']*l2**2), ('K3m*l2*lf', nv['K3m']*l2*lf), ('K3m*lf^2', nv['K3m']*lf**2)]
B += [('L', nv['Lchi5_5']), ('r5*L', r5*nv['Lchi5_5'])]
log(f"step_05: {len(B)} base atoms ({time.time()-t0:.0f}s)")

def find_relation(allvals, dps):
    mp.mp.dps = dps + 30
    n = len(allvals)
    C = mp.mpf(10) ** (dps - 30)
    M = fmpz_mat(n, n + 1)
    for i, v in enumerate(allvals):
        M[i, i] = fmpz(1)
        M[i, n] = fmpz(int(mp.nint(C * v)))
    R = M.lll(delta=0.99, eta=0.51, gram='exact')
    mp.mp.dps = DPS + 40
    return [[int(R[r, cc]) for cc in range(n + 1)] for r in range(n)]

def normalize(row):
    coeffs = [int(cc) for cc in row[:-1]]
    g = reduce(gcd, [abs(cc) for cc in coeffs if cc]) or 1
    coeffs = [cc // g for cc in coeffs]
    if next((cc for cc in coeffs if cc), 1) < 0: coeffs = [-cc for cc in coeffs]
    return coeffs

# ---- planted-relation control: target Li5(rho^2) with duplication relation present
cols_ctl = [nv['Li5(rho)'], nv['Li5(-rho)'], nv['Li5(rho^6)'], l2, lf, z5]
rows = find_relation([nv['Li5(rho^2)'] + cols_ctl, 300)
planted = None
for row in rows:
    if row[0] == 16 and row[1] == 16 and row[2] == 16 and row[3] == -1:
        planted = row
log(f"planted control (16Li5(rho)+16Li5(-rho)-Li5(rho^2)=0): {'FOUND' if planted else 'NOT FOUND (setup issue)'}"
    if planted else f"planted control: rows[0]={rows[0][:6]}")
assert planted

names = [n for n, _ in B] + ['r5*' + n for n, _ in B]
cols = [v for _, v in B] + [r5*v for _, v in B]
result = {}
for tag in ('gol1', 'gol3'):
    target = nv['R_' + tag]
    cands = {}
    for dps in (1150, 1400):
        rows = find_relation([target] + cols, dps)
        best = None
        for row in rows:
            if row[0] != 0:
                h = max(abs(x) for x in row[:-1])
                if best is None or h < best[0]: best = (h, normalize(row))
        cands[dps] = best
        h, row = best
        res = Fr(row[0]) * target
        for c, v in zip(row[1:], cols):
            res += Fr(c) * v
        log(f"  {tag} dps {dps}: height {h}, residual@1490 {mp.nstr(abs(res), 3)}")
    same = cands[1150][1] == cands[1400][1]
    log(f"  {tag}: identical vectors at 1150/1400: {same}")
    row = cands[1400][1]
    Li = names.index('L'); rLi = names.index('r5*L')
    log(f"  {tag}: target coeff {row[0]}; L coeff {row[1+Li]}; sqrt5*L coeff {row[1+rLi]}")
    log(f"     => coefficient of L(chi5,5) in R: ({-row[1+Li]}/{row[0]}) + ({row[1+rLi]}/{row[0]})*sqrt5")
    result[tag] = {'relation': row, 'names': ['TARGET'] + names,
                   'same_two_prec': bool(same), 'height': max(abs(x) for x in row),
                   'L_coeff': [-row[1+Li], row[0]], 'r5L_coeff': [row[1+rLi], row[0]]}
    # residual at full precision for the record
    res = Fr(row[0]) * target
    for c, v in zip(row[1:], cols):
        res += Fr(c) * v
    result[tag]['residual_1490'] = mp.nstr(abs(res), 3)
json.dump(result, open(GOLD + 'remainder_fits.json', 'w'), indent=1)
log(f"step_05 fits done ({time.time()-t0:.0f}s); saved remainder_fits.json")
'''
open('step_05.py', 'w').write(step05)
import subprocess
r = subprocess.run(['nice', '-n', '10', '<project>/.venv/bin/python', 'step_05.py'], capture_output=True, text=True, cwd='.', timeout=2400)
print(r.stdout)
print("STDERR:", r.stderr[-1500:])
