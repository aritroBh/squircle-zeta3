import json, mpmath as mp
from fractions import Fraction as Fr
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
DPS = 1490
mp.mp.dps = DPS + 40
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
names = [n for n, _ in B]
vals = [v for _, v in B]
print(f"{len(B)} atoms; with sqrt5 columns: {2*len(B)}")
from flint import fmpz_mat, fmpz
def find_relation(allvals, dps, C_exp=None):
    mp.mp.dps = dps + 30
    n = len(allvals)
    C = mp.mpf(10) ** (dps - 30 if C_exp is None else C_exp)
    M = fmpz_mat(n, n + 1)
    for i, v in enumerate(allvals):
        M[i, i] = fmpz(1)
        M[i, n] = fmpz(int(mp.nint(C * v)))
    R = M.lll(delta=0.99, eta=0.51, gram='exact')
    mp.mp.dps = DPS + 40
    return [[int(R[r, cc]) for cc in range(n + 1)] for r in range(n)]
# quick probe at 200 digits: target R_gol1, columns = atoms + sqrt5*atoms
cols = []
colnames = []
for nm, v in B:
    cols.append(v); colnames.append(nm)
for nm, v in B:
    cols.append(r5*v); colnames.append('r5*' + nm)
probe = [nv['R_gol1']] + cols
rows = find_relation(probe, 200)
for row in rows[:4]:
    h = max(abs(x) for x in row[:-1])
    resid_scaled = row[-1]
    if row[0] != 0 or h < 10**6:
        print("row0 h=%.3g residcol=%.3g:" % (h, abs(resid_scaled)), row[0], row[1:8], "...")
