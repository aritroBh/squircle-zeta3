
step09 = r'''"""step_09.py -- G2b: classical closed form, reduced by the Lie_5 class.

class(R4^(5)(1/2) - N5) = (359/2, 0, 275/2, 0, 247/2, 122) lies in the Li5 span (rank-4 check,
exact). A particular Li5 representative (kernel taus set to 0):
  P5 = -14*Li5(-2) - (11/2)*Li5(1/3) + (219/2)*Li5(2/3) + (7/2)*Li5(3/4).
Search C: T2 = R4^(5)(1/2) - N5 - P5 over [weight-5 product basis + zeta(5)] at two precisions.
Uses the dps-2100 values (same code as step_08).
"""
import json, re, time, mpmath as mp
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
from flint import fmpz_mat, fmpz

BASE = "<project>/research/08_theory"
G2B = BASE + "/G2b_point_depth_M4/"
G2 = BASE + "/G2_M4alpha/"
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
log(f"R4^(5)(1/2) dps {DPS} ({time.time()-t0:.0f}s)")

def li32(x, y):
    x = mp.mpf(x); y = mp.mpf(y)
    Nl = int((DPS + 40) * mp.log(10) / -mp.log(abs(x*y))) + 100
    S = mp.mpf(0); ym = y; xn = x; tot = mp.mpf(0)
    for n in range(2, Nl + 1):
        S += ym / ((n-1)*(n-1)); ym *= y; xn *= x
        tot += xn * S / (n*n*n)
    return tot
A1 = li32(mp.mpf(1)/4, mp.mpf(4)/3)
A2 = li32(mp.mpf(8)/9, mp.mpf(1))
log(f"atoms done ({time.time()-t0:.0f}s)")

def Li(s, x): return mp.polylog(s, mp.mpf(x))
l2, l3 = mp.log(2), mp.log(3)
pi2, pi4 = mp.pi**2, mp.pi**4
z3, z5 = mp.zeta(3), mp.zeta(5)
Li2a = Li(2, mp.mpf(1)/3)
L313, L31m = Li(3, mp.mpf(1)/3), Li(3, -mp.mpf(1)/2)
L413, L423, L412, L419, L434 = Li(4, mp.mpf(1)/3), Li(4, mp.mpf(2)/3), Li(4, mp.mpf(1)/2), Li(4, mp.mpf(1)/9), Li(4, mp.mpf(3)/4)

P5 = (-14)*Li(5, -2) + Fr(-11, 2)*Li(5, mp.mpf(1)/3) + Fr(219, 2)*Li(5, mp.mpf(2)/3) + Fr(7, 2)*Li(5, mp.mpf(3)/4)
N5V = Fr(5, 4)*A1 + Fr(1, 3)*A2
T2 = TARGET - N5V - P5
log(f"T2 = R4^(5) - N5 - P5 assembled ({time.time()-t0:.0f}s)")

# product basis (weight-5 products + z5), no Li5 primitives
BASIS = [('z5', z5), ('z3*Li2(1/3)', z3*Li2a), ('z3*pi2', z3*pi2)]
for nm, vv in (('l2^2', l2**2), ('l2*l3', l2*l3), ('l3^2', l3**2)):
    BASIS.append((f'z3*{nm}', z3*vv))
BASIS += [('pi2*Li3(1/3)', pi2*L313), ('pi2*Li3(-1/2)', pi2*L31m)]
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
    BASIS.append((f'({nm})*l2', vv*l2))
    BASIS.append((f'({nm})*l3', vv*l3))
w3basis = [z3, L313, L31m, Li2a*l2, Li2a*l3, pi2*l2, pi2*l3, l2**3, l2**2*l3, l2*l3**2, l3**3]
w3names = ['z3', 'Li3(1/3)', 'Li3(-1/2)', 'Li2*l2', 'Li2*l3', 'pi2*l2', 'pi2*l3',
           'l2^3', 'l2^2l3', 'l2l3^2', 'l3^3']
for nm, vv in zip(w3names, w3basis):
    BASIS.append((f'({nm})*l2^2', vv*l2**2))
    BASIS.append((f'({nm})*l2*l3', vv*l2*l3))
    BASIS.append((f'({nm})*l3^2', vv*l3**2))
for nm, vv in (('Li2', Li2a), ('pi2', pi2)):
    for ln, lv in (('l2^3', l2**3), ('l2^2l3', l2**2*l3), ('l2l3^2', l2*l3**2), ('l3^3', l3**3)):
        BASIS.append((f'{nm}*{ln}', vv*lv))
for ln, lv in (('l2^4', l2**4), ('l2^3l3', l2**3*l3), ('l2^2l3^2', l2**2*l3**2),
               ('l2l3^3', l2*l3**3), ('l3^4', l3**4)):
    BASIS.append((f'l2*{ln}', l2*lv))
    BASIS.append((f'l3*{ln}', l3*lv))
for i in range(6):
    BASIS.append((f'l2^{5-i}*l3^{i}', l2**(5-i)*l3**i))
BASIS.append(('Li3(1/3)*Li2(1/3)', L313*Li2a))
BASIS.append(('Li3(-1/2)*Li2(1/3)', L31m*Li2a))
# kernel Li5 combinations (zero class; their values are products + z5)
K5 = [('Li5(1/2)', Li(5, mp.mpf(1)/2)),
      ('Li5(1/9)-48*Li5(1/3)', Li(5, mp.mpf(1)/9) - 48*Li(5, mp.mpf(1)/3)),
      ('Li5(-1/3)-2*Li5(1/3)', Li(5, -mp.mpf(1)/3) - 2*Li(5, mp.mpf(1)/3)),
      ('Li5(-1/2)-Li5(-2)', Li(5, -mp.mpf(1)/2) - Li(5, -2)),
      ('Li5(1/4)-16*Li5(-2)', Li(5, mp.mpf(1)/4) - 16*Li(5, -2)),
      ('Li5(-1/8)-162*Li5(-2)', Li(5, -mp.mpf(1)/8) - 162*Li(5, -2))]
BASIS += K5
seen = set(); ded = []
for nm, vv in BASIS:
    if nm in seen: continue
    seen.add(nm); ded.append((nm, vv))
BASIS = ded
log(f"product basis: {len(BASIS)} elements ({time.time()-t0:.0f}s)")

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
    return [[R[r, cc] for cc in range(n + 1)] for r in range(n)]

def normalize(row):
    coeffs = [int(cc) for cc in row[:-1]]
    g = reduce(gcd, [abs(cc) for cc in coeffs if cc]) or 1
    coeffs = [cc // g for cc in coeffs]
    if next((cc for cc in coeffs if cc), 1) < 0: coeffs = [-cc for cc in coeffs]
    return coeffs

names = [nm for nm, _ in BASIS]
vals = [v for _, v in BASIS]
cands = {}
for dps in (1500, 1900):
    rows = find_relation([T2] + vals, dps)
    row = next((r for r in rows if r[0] != 0), None)
    cands[dps] = normalize(row)
    resid = sum(cc * v for cc, v in zip(cands[dps], [T2] + vals))
    h = max(abs(cc) for cc in cands[dps])
    log(f"  C dps {dps}: height {h}, resid@full = {mp.nstr(abs(resid), 3)}")
same = cands[1500] == cands[1900]
log(f"  C: identical at two precisions: {same}")
if same:
    a = cands[1900]
    terms_ = [(Fr(-cc, a[0]), nm) for cc, nm in zip(a[1:], names) if cc]
    log(f"  C relation ({len(terms_)} terms, a0 = {a[0]}):")
    for cf, nm in terms_:
        log(f"      {cf!s:>12} * {nm}")
    json.dump(dict(a0=a[0], coeffs=a[1:], names=names), open(G2B + 'closed_C.json', 'w'), indent=1)
log(f"step_09 done ({time.time()-t0:.0f}s)")
'''
G2B = "<project>/research/08_theory/G2b_point_depth_M4/"
open(G2B + "step_09.py", "w").write(step09)
import subprocess, time
PY = "<project>/.venv/bin/python"
t0 = time.time()
r = subprocess.run(["nice", "-n", "19", PY, "step_09.py"], cwd=G2B, capture_output=True, text=True, timeout=1800)
print("elapsed", time.time() - t0, "rc", r.returncode)
print(r.stdout)
print("STDERR:", r.stderr[-2500:])
