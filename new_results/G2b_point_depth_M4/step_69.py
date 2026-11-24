
step06 = r'''"""step_06.py -- G2b Task 4 (numerical): closed form of R4^{(5)}(1/2) over MT(Z[1/6]) at weight 5
plus the cobracket-fixed depth-two correction N5 = (3/4)*Li32(8/9,3/4) + 4*Li32(2,1/3).

Two searches, each at two lattice precisions, values from values_hex.json (1300 digits, step_04):
  (A) target R4^{(5)}(1/2) over [classical basis] + [the two Li32 atoms, free coefficients]:
      the LLL must recover coefficients 3/4 and 4 (independent confirmation of the cobracket).
  (B) target R4^{(5)}(1/2) - N5 over [classical basis] only.
Residuals printed at full 1300 digits.
"""
import json, time, mpmath as mp
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
from flint import fmpz_mat, fmpz

BASE = "<project>/research/08_theory"
G2B = BASE + "/G2b_point_depth_M4/"
DPS = 1300
mp.mp.dps = DPS + 30
t0 = time.time()

def log(msg):
    print(msg, flush=True)
    with open(G2B + "progress.log", "a") as f:
        f.write(msg + "\n")

V = {k: mp.mpf(v) for k, v in json.load(open(G2B + "values_hex.json")).items()}
l2, l3 = mp.log(2), mp.log(3)
pi2, pi4 = mp.pi**2, mp.pi**4
z3, z5 = mp.zeta(3), mp.zeta(5)
Li2a = V['Li2(1/3)']
L313, L31m = V['Li3(1/3)'], V['Li3(-1/2)']
L413, L423, L412, L419, L434 = V['Li4(1/3)'], V['Li4(2/3)'], V['Li4(1/2)'], V['Li4(1/9)'], V['Li4(3/4)']

# ---- classical weight-5 basis (COR5-style reduced spanning set)
BASIS = []
def add(name, val):
    BASIS.append((name, val))

for a in ['-2', '1/3', '2/3', '1/2', '1/9', '-1/2', '-1/3', '1/4', '3/4', '-1/8', '8/9', '-1']:
    add(f'Li5({a})', V[f'Li5({a})'])
add('z5', z5)
add('z3*Li2(1/3)', z3*Li2a)
add('z3*pi2', z3*pi2)
for nm, vv in (('l2^2', l2**2), ('l2*l3', l2*l3), ('l3^2', l3**2)):
    add(f'z3*{nm}', z3*vv)
add('pi2*Li3(1/3)', pi2*L313)
add('pi2*Li3(-1/2)', pi2*L31m)
# COR5 25-element weight-4 basis x logs (subset: the non-product Li4's and key products)
w4basis = [L413, L423, L434, L412, L419,
           L313*l2, L313*l3, L31m*l2, L31m*l3, z3*l2, z3*l3,
           Li2a**2, Li2a*l2**2, Li2a*l2*l3, Li2a*l3**2, pi2*Li2a,
           pi2*l2**2, pi2*l2*l3, pi2*l3**2, pi4,
           l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4]
w4names = ['Li4(1/3)', 'Li4(2/3)', 'Li4(3/4)', 'Li4(1/2)', 'Li4(1/9)',
           'Li3(1/3)l2', 'Li3(1/3)l3', 'Li3(-1/2)l2', 'Li3(-1/2)l3', 'z3l2', 'z3l3',
           'Li2^2', 'Li2*l2^2', 'Li2*l2l3', 'Li2*l3^2', 'pi2*Li2',
           'pi2*l2^2', 'pi2*l2l3', 'pi2*l3^2', 'pi4',
           'l2^4', 'l2^3l3', 'l2^2l3^2', 'l2l3^3', 'l3^4']
for nm, vv in zip(w4names, w4basis):
    add(f'({nm})*l2', vv*l2)
    add(f'({nm})*l3', vv*l3)
# weight-3 basis x log^2 monomials
w3basis = [z3, L313, L31m, Li2a*l2, Li2a*l3, pi2*l2, pi2*l3, l2**3, l2**2*l3, l2*l3**2, l3**3]
w3names = ['z3', 'Li3(1/3)', 'Li3(-1/2)', 'Li2*l2', 'Li2*l3', 'pi2*l2', 'pi2*l3',
           'l2^3', 'l2^2l3', 'l2l3^2', 'l3^3']
for nm, vv in zip(w3names, w3basis):
    add(f'({nm})*l2^2', vv*l2**2)
    add(f'({nm})*l2*l3', vv*l2*l3)
    add(f'({nm})*l3^2', vv*l3**2)
# weight-2 x log^3, logs x log^4, log^5
for nm, vv in (('Li2', Li2a), ('pi2', pi2)):
    for ln, lv in (('l2^3', l2**3), ('l2^2l3', l2**2*l3), ('l2l3^2', l2*l3**2), ('l3^3', l3**3)):
        add(f'{nm}*{ln}', vv*lv)
for ln, lv in (('l2^4', l2**4), ('l2^3l3', l2**3*l3), ('l2^2l3^2', l2**2*l3**2),
               ('l2l3^3', l2*l3**3), ('l3^4', l3**4)):
    add(f'l2*{ln}', l2*lv)
    add(f'l3*{ln}', l3*lv)
for i in range(6):
    add(f'l2^{5-i}*l3^{i}', l2**(5-i)*l3**i)
add('Li3(1/3)*Li2(1/3)', L313*Li2a)
add('Li3(-1/2)*Li2(1/3)', L31m*Li2a)
add('Li2(1/3)^2*l2', Li2a**2*l2)
add('Li2(1/3)^2*l3', Li2a**2*l3)
log(f"classical basis: {len(BASIS)} elements ({time.time()-t0:.0f}s)")

ATOMS = [('Li32(8/9,3/4)', V['Li32(8/9,3/4)']), ('Li32(2,1/3)', V['Li32(2,1/3)'])]
TARGET = V['R4_w5']
N5VAL = Fr(3, 4) * 0  # placeholder

def find_relation(vals, dps, margin=40, delta=0.99):
    mp.mp.dps = dps
    C = mp.mpf(10) ** (dps - margin)
    n = len(vals)
    M = fmpz_mat(n, n + 1)
    for i, v in enumerate(vals):
        M[i, i] = fmpz(1)
        M[i, n] = fmpz(int(mp.nint(C * v)))
    R = M.lll(delta=delta, eta=0.51)
    mp.mp.dps = DPS + 30
    return [[R[r, c] for c in range(n + 1)] for r in range(n)]

def normalize(row):
    coeffs = [int(c) for c in row[:-1]]
    g = reduce(gcd, [abs(c) for c in coeffs if c]) or 1
    coeffs = [c // g for c in coeffs]
    if next((c for c in coeffs if c), 1) < 0: coeffs = [-c for c in coeffs]
    return coeffs

def search(tag, target, items, precisions):
    names = [nm for nm, _ in items]
    vals = [v for _, v in items]
    cands = {}
    for dps in precisions:
        rows = find_relation([target] + vals, dps)
        cands[dps] = normalize(rows[0])
        resid = sum(c * v for c, v in zip(cands[dps], [target] + vals))
        h = max(abs(c) for c in cands[dps])
        log(f"  {tag} dps {dps}: height {h}, resid@1300 = {mp.nstr(abs(resid), 3)}")
    same = cands[precisions[0]] == cands[precisions[1]]
    log(f"  {tag}: identical at two precisions: {same}")
    if same:
        a = cands[precisions[1]]
        terms = [(Fr(-c, a[0]), nm) for c, nm in zip(a[1:], names) if c]
        log(f"  {tag} relation ({len(terms)} terms, a0 = {a[0]}):")
        for cf, nm in terms:
            log(f"      {cf!s:>12} * {nm}")
        return dict(a0=a[0], coeffs=a[1:], names=names, target=tag)
    return None

if __name__ == '__main__':
    log(f"step_06 start {time.strftime('%H:%M:%S')}, dps {DPS}")
    # (A) free Li32 coefficients: must come out 3/4 and 4
    resA = search('A: R4^(5)(1/2) over classical + 2 Li32 atoms', TARGET, BASIS + ATOMS, (800, 1100))
    if resA:
        i1 = resA['names'].index('Li32(8/9,3/4)'); i2 = resA['names'].index('Li32(2,1/3)')
        c1 = Fr(-resA['coeffs'][i1], resA['a0']); c2 = Fr(-resA['coeffs'][i2], resA['a0'])
        log(f"  A: Li32 coefficients from LLL: {c1}, {c2} (cobracket-fixed: 3/4, 4): "
            f"{'MATCH' if (c1, c2) == (Fr(3, 4), Fr(4)) else 'MISMATCH'}")
        resA['li32_coeffs'] = [str(c1), str(c2)]
        json.dump(resA, open(G2B + 'closed_A.json', 'w'), indent=1)
    # (B) N5 fixed: R4^(5) - N5 over classical only
    n5 = Fr(3, 4) * ATOMS[0][1] + Fr(4) * ATOMS[1][1]
    resB = search('B: R4^(5)(1/2) - N5 over classical', TARGET - n5, BASIS, (800, 1100))
    if resB:
        json.dump(resB, open(G2B + 'closed_B.json', 'w'), indent=1)
    log(f"step_06 done ({time.time()-t0:.0f}s)")
'''
G2B = "<project>/research/08_theory/G2b_point_depth_M4/"
open(G2B + "step_06.py", "w").write(step06)
import subprocess, time
PY = "<project>/.venv/bin/python"
t0 = time.time()
r = subprocess.run(["nice", "-n", "19", PY, "step_06.py"], cwd=G2B, capture_output=True, text=True, timeout=1800)
print("elapsed", time.time() - t0, "rc", r.returncode)
print(r.stdout)
print("STDERR:", r.stderr[-2500:])
