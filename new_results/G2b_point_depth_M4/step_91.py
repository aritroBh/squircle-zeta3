
step08 = r'''"""step_08.py -- G2b: classical closed form of R4^(5)(1/2) and of R4^(5)(1/2) - N5, at dps 2100.

Recomputes the needed values at dps 2100 (HPL series at c=1/2; Li32 single series) and runs the
two LLL searches of step_06 at lattice precisions (1700, 2000). Detection threshold with
n = 133 numbers is log10(height) ~ dps/132; dps 1300 was below the candidate heights.
"""
import json, re, time, mpmath as mp
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
from flint import fmpz_mat, fmpz

BASE = "<project>/research/08_theory"
G2B = BASE + "/G2b_point_depth_M4/"
G2 = BASE + "/G2_M4alpha/"
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
DPS = 2100
mp.mp.dps = DPS + 40
t0 = time.time()

def log(msg):
    print(msg, flush=True)
    with open(G2B + "progress.log", "a") as f:
        f.write(msg + "\n")

def hpl(w, x, N):
    a = w[-1]
    s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1):
        s[m] = mp.mpf(1) / m if a == 1 else mp.mpf((-1) ** (m + 1)) / m
    for a in reversed(w[:-1]):
        t = [mp.mpf(0)] * (N + 1)
        if a == 0:
            for m in range(1, N + 1): t[m] = s[m] / m
        else:
            acc = mp.mpf(0)
            for m in range(0, N):
                acc = acc * a + s[m]; t[m + 1] = acc / (m + 1)
        s = t
    return mp.fsum(s[m] * x ** m for m in range(1, N + 1))

def terms(f):
    out = []
    for nm, v in json.load(open(f)):
        q = Fr(v)
        if nm == '1': out.append(('1', (), 0, q)); continue
        m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]', nm)
        w = tuple(int(t) for t in m.group(3).split(',')) if m.group(3) else ()
        out.append((m.group(2), w, int(m.group(1)), q))
    return out

R4 = terms(G2 + "hpl_fit_k4_(0,0,0).json")
c = mp.mpf(1) / 2
N = int((DPS + 25) / mp.log10(2)) + 60
log(f"step_08: dps {DPS}, series N = {N} ({time.time()-t0:.0f}s)")
cache = {}
def H(w, x):
    if (w, x) not in cache: cache[(w, x)] = hpl(w, x, N)
    return cache[(w, x)]
TARGET = mp.mpf(0)
for kind, w, e, q in R4:
    if kind == '1' or len(w) != 5: continue
    qf = mp.mpf(q.numerator) / q.denominator
    f = c ** e * H(w, c)
    if kind == 'He': f += (-c) ** e * H(w, -c)
    TARGET += qf * f
log(f"R4^(5)(1/2) at dps {DPS} done ({time.time()-t0:.0f}s)")

def li32(x, y):
    x = mp.mpf(x); y = mp.mpf(y)
    r = abs(x * y)
    Nl = int((DPS + 40) * mp.log(10) / -mp.log(r)) + 100
    S = mp.mpf(0); ym = y; xn = x; tot = mp.mpf(0)
    for n in range(2, Nl + 1):
        S += ym / ((n-1)*(n-1)); ym *= y; xn *= x
        tot += xn * S / (n*n*n)
    return tot
ATOMS = [('Li32(1/4,4/3)', li32(mp.mpf(1)/4, mp.mpf(4)/3)), ('Li32(8/9,1)', li32(mp.mpf(8)/9, mp.mpf(1)))]
log(f"atoms done ({time.time()-t0:.0f}s)")

def Li(s, x): return mp.polylog(s, mp.mpf(x))
l2, l3 = mp.log(2), mp.log(3)
pi2, pi4 = mp.pi**2, mp.pi**4
z3, z5 = mp.zeta(3), mp.zeta(5)
Li2a = Li(2, mp.mpf(1)/3)
L313, L31m = Li(3, mp.mpf(1)/3), Li(3, -mp.mpf(1)/2)
L413, L423, L412, L419, L434 = Li(4, mp.mpf(1)/3), Li(4, mp.mpf(2)/3), Li(4, mp.mpf(1)/2), Li(4, mp.mpf(1)/9), Li(4, mp.mpf(3)/4)

BASIS = []
def add(name, val): BASIS.append((name, val))
for a in ['-2', '1/3', '2/3', '1/2', '1/9', '-1/2', '-1/3', '1/4', '3/4', '-1/8', '8/9']:
    x = eval(a, {}, {}) if False else Fr(a)
    add(f'Li5({a})', Li(5, x.numerator/x.denominator))
add('z5', z5)
add('z3*Li2(1/3)', z3*Li2a)
add('z3*pi2', z3*pi2)
for nm, vv in (('l2^2', l2**2), ('l2*l3', l2*l3), ('l3^2', l3**2)):
    add(f'z3*{nm}', z3*vv)
add('pi2*Li3(1/3)', pi2*L313)
add('pi2*Li3(-1/2)', pi2*L31m)
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
w3basis = [z3, L313, L31m, Li2a*l2, Li2a*l3, pi2*l2, pi2*l3, l2**3, l2**2*l3, l2*l3**2, l3**3]
w3names = ['z3', 'Li3(1/3)', 'Li3(-1/2)', 'Li2*l2', 'Li2*l3', 'pi2*l2', 'pi2*l3',
           'l2^3', 'l2^2l3', 'l2l3^2', 'l3^3']
for nm, vv in zip(w3names, w3basis):
    add(f'({nm})*l2^2', vv*l2**2)
    add(f'({nm})*l2*l3', vv*l2*l3)
    add(f'({nm})*l3^2', vv*l3**2)
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
seen = set(); ded = []
for nm, vv in BASIS:
    if nm in seen: continue
    seen.add(nm); ded.append((nm, vv))
BASIS = ded
log(f"basis: {len(BASIS)} + 2 atoms ({time.time()-t0:.0f}s)")

def find_relation(vals, dps, margin=40, delta=0.99):
    mp.mp.dps = dps
    C = mp.mpf(10) ** (dps - margin)
    n = len(vals)
    M = fmpz_mat(n, n + 1)
    for i, v in enumerate(vals):
        M[i, i] = fmpz(1)
        M[i, n] = fmpz(int(mp.nint(C * v)))
    R = M.lll(delta=delta, eta=0.51)
    mp.mp.dps = DPS + 40
    return [[R[r, c2] for c2 in range(n + 1)] for r in range(n)]

def normalize(row):
    coeffs = [int(c2) for c2 in row[:-1]]
    g = reduce(gcd, [abs(c2) for c2 in coeffs if c2]) or 1
    coeffs = [c2 // g for c2 in coeffs]
    if next((c2 for c2 in coeffs if c2), 1) < 0: coeffs = [-c2 for c2 in coeffs]
    return coeffs

def search(tag, target, items, precisions):
    names = [nm for nm, _ in items]
    vals = [v for _, v in items]
    cands = {}
    for dps in precisions:
        rows = find_relation([target] + vals, dps)
        row = next((r for r in rows if r[0] != 0), None)
        if row is None:
            log(f"  {tag} dps {dps}: no relation involving the target"); return None
        cands[dps] = normalize(row)
        resid = sum(c2 * v for c2, v in zip(cands[dps], [target] + vals))
        h = max(abs(c2) for c2 in cands[dps])
        log(f"  {tag} dps {dps}: height {h}, resid@full = {mp.nstr(abs(resid), 3)}")
    same = cands[precisions[0]] == cands[precisions[1]]
    log(f"  {tag}: identical at two precisions: {same}")
    if same:
        a = cands[precisions[1]]
        terms_ = [(Fr(-c2, a[0]), nm) for c2, nm in zip(a[1:], names) if c2]
        log(f"  {tag} relation ({len(terms_)} terms, a0 = {a[0]}):")
        for cf, nm in terms_:
            log(f"      {cf!s:>12} * {nm}")
        return dict(a0=a[0], coeffs=a[1:], names=names, target=tag)
    return None

if __name__ == '__main__':
    log(f"step_08 search start ({time.time()-t0:.0f}s)")
    resA = search('A: R4^(5)(1/2) over classical + 2 Li32 atoms', TARGET, BASIS + ATOMS, (1700, 2000))
    if resA:
        i1 = resA['names'].index('Li32(1/4,4/3)'); i2 = resA['names'].index('Li32(8/9,1)')
        c1 = Fr(-resA['coeffs'][i1], resA['a0']); c2 = Fr(-resA['coeffs'][i2], resA['a0'])
        log(f"  A: Li32 coefficients from LLL: {c1}, {c2} (target 5/4, 1/3): "
            f"{'MATCH' if (c1, c2) == (Fr(5, 4), Fr(1, 3)) else 'MISMATCH'}")
        resA['li32_coeffs'] = [str(c1), str(c2)]
        json.dump(resA, open(G2B + 'closed_A.json', 'w'), indent=1)
    n5 = Fr(5, 4) * ATOMS[0][1] + Fr(1, 3) * ATOMS[1][1]
    resB = search('B: R4^(5)(1/2) - N5 over classical', TARGET - n5, BASIS, (1700, 2000))
    if resB:
        json.dump(resB, open(G2B + 'closed_B.json', 'w'), indent=1)
    log(f"step_08 done ({time.time()-t0:.0f}s)")
'''
G2B = "<project>/research/08_theory/G2b_point_depth_M4/"
open(G2B + "step_08.py", "w").write(step08)
import subprocess, time
PY = "<project>/.venv/bin/python"
t0 = time.time()
r = subprocess.run(["nice", "-n", "19", PY, "step_08.py"], cwd=G2B, capture_output=True, text=True, timeout=1800)
print("elapsed", time.time() - t0, "rc", r.returncode)
print(r.stdout)
print("STDERR:", r.stderr[-2500:])
