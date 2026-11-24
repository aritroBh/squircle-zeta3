import json, mpmath as mp
from fractions import Fraction as Fr
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
DPS = 1490
mp.mp.dps = DPS + 40
nv = {k: mp.mpf(v) for k, v in json.load(open(GOLD + 'number_values.json')).items()}
r5 = mp.sqrt(5)
phi = (1 + r5)/2
rho = 1/phi
# add Li4 primitives and Li2 primitives
add = {}
def Li(s, x): return mp.polylog(s, mp.mpf(x))
add['Li4(-1/(2phi))'] = Li(4, -1/(2*phi))
add['Li4(1/(2phi^2))'] = Li(4, 1/(2*phi**2))
add['Li4(2phi^-2)'] = Li(4, 2*phi**-2)
add['Li4(phi/2)'] = Li(4, phi/2)
add['Li4(-rho^3)'] = Li(4, -rho**3)
add['Li4(rho^3)'] = Li(4, rho**3)
add['Li2(rho^3)'] = Li(2, rho**3)
nv.update(add)
json.dump({k: mp.nstr(v, DPS) for k, v in nv.items()}, open(GOLD + 'number_values.json', 'w'), indent=1)
print("added Li4/Li2 primitives")

# build enlarged basis and probe at 600 digits
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
# invisible golden constants
B += [('Li5(-rho)', nv['Li5(-rho)']), ('Li5(rho^2)', nv['Li5(rho^2)'])]
B += [('K4m*l2', nv['K4m']*l2), ('K4m*lf', nv['K4m']*lf), ('K4q*l2', nv['K4q']*l2), ('K4q*lf', nv['K4q']*lf)]
B += [('K3m*l2^2', nv['K3m']*l2**2), ('K3m*l2*lf', nv['K3m']*l2*lf), ('K3m*lf^2', nv['K3m']*lf**2)]
# Li4 primitives x logs
for nm in ['Li4(-1/(2phi))', 'Li4(1/(2phi^2))', 'Li4(2phi^-2)', 'Li4(phi/2)', 'Li4(-rho^3)', 'Li4(rho^3)']:
    B.append((nm + '*l2', nv[nm]*l2))
    B.append((nm + '*lf', nv[nm]*lf))
# Li2 primitive x weight-3 and powers
L2b = nv['Li2(rho^3)']
B += [('Li2(rho^3)*z3', L2b*z3), ('Li2(rho^3)*K3m', L2b*nv['K3m'])]
for nm_, v_ in [('l2^3', l2**3), ('l2^2*lf', l2**2*lf), ('l2*lf^2', l2*lf**2), ('lf^3', lf**3)]:
    B.append((f'Li2(rho^3)*{nm_}', L2b*v_))
B += [('Li2(rho^3)*pi2*l2', L2b*pi2*l2), ('Li2(rho^3)*pi2*lf', L2b*pi2*lf)]
B += [('Li2(rho^3)^2*l2', L2b**2*l2), ('Li2(rho^3)^2*lf', L2b**2*lf)]
B += [('L', nv['Lchi5_5'])]
print(f"{len(B)} base atoms")
from flint import fmpz_mat, fmpz
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
cols = [v for _, v in B] + [r5*v for _, v in B]
import time
for tag in ('gol1', 'gol3'):
    t0 = time.time()
    rows = find_relation([nv['R_' + tag]] + cols, 600)
    best = None
    for row in rows:
        if row[0] != 0:
            h = max(abs(x) for x in row[:-1])
            if best is None or h < best[0]: best = (h, row)
    if best:
        h, row = best
        res = Fr(row[0]) * nv['R_' + tag]
        for c, v in zip(row[1:], cols):
            res += Fr(c) * v
        print(f"{tag} @600d: height {h}, residual@1490 {mp.nstr(abs(res), 3)} ({time.time()-t0:.0f}s)")
    else:
        print(f"{tag} @600d: none ({time.time()-t0:.0f}s)")
