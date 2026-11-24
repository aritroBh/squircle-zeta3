import json, mpmath as mp
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
cols = []
for nm, v in B: cols.append(v)
for nm, v in B: cols.append(r5*v)
for tag in ('gol1', 'gol3'):
    probe = [nv['R_' + tag]] + cols
    rows = find_relation(probe, 400)
    best = None
    for row in rows:
        if row[0] != 0:
            h = max(abs(x) for x in row[:-1])
            if best is None or h < best[0]: best = (h, row)
    if best:
        resid = sum(Fr(c)*mp.mpf(0) for c in []) # placeholder
        print(f"{tag} @400d: target relation found, height {best[0]}")
        # evaluate residual at full precision
        res = Fr(best[0][0]) * nv['R_' + tag]
        for c, v in zip(best[0][1:], cols):
            res += Fr(c) * v
        print(f"   residual at 1490d: {mp.nstr(abs(res), 3)}")
    else:
        print(f"{tag} @400d: no target relation in first rows")
