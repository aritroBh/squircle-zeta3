"""step_04.py -- G2b numerics: R4's weight-5 part and M4 at the hexagon corner (alpha = 2pi/3,
c = 1/2) to ~1280 digits, plus Li_{3,2} atoms at {2,3}-units and Li3/Li5 reference values.

Method: planner_check_b4.py's independent mpmath HPL series evaluator (Remiddi-Vermaseren
recursion); at c = 1/2 the series converges like 2^-N. Li_{3,2}(x,y) = sum_{n>m>0} x^n y^m/(n^3 m^2)
as a single series with incremental inner sum S_y(n) = sum_{m<n} y^m/m^2 (exact at each n).
"""
import json, re, time, mpmath as mp
from fractions import Fraction as Fr

BASE = "<project>/research/08_theory"
G2B = BASE + "/G2b_point_depth_M4/"
G2 = BASE + "/G2_M4alpha/"
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
DPS = 1300
mp.mp.dps = DPS + 30
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
P2 = terms(G2 + "hpl_fit_k4_(1,0,0).json")
PS = terms(G1 + "hpl_fit_k3_(1,0,0).json")

c = mp.mpf(1) / 2
N = int((DPS + 25) / mp.log10(2)) + 50
log(f"step_04: dps {DPS}, series order N = {N} ({time.time()-t0:.0f}s)")

cache = {}
def H(w, x):
    key = (w, x)
    if key not in cache:
        cache[key] = hpl(w, x, N)
    return cache[key]

def ev(ts, weight=None):
    tot = mp.mpf(0)
    for kind, w, e, q in ts:
        if kind == '1':
            if weight in (None, 0): tot += mp.mpf(q.numerator) / q.denominator
            continue
        if weight is not None and len(w) != weight: continue
        qf = mp.mpf(q.numerator) / q.denominator
        f = c ** e * H(w, c)
        if kind == 'He':
            # He[w] = c^e H_w(c) + (-c)^e H_w(-c): evaluate H_w at -c (planner_check_b4).
            # (Equivalently H_w(-c) = (-1)^nz(w) H_{-w}(c); the earlier version missed (-1)^nz.)
            f += (-c) ** e * H(w, -c)
        tot += qf * f
    return tot

vals = {}
vals['R4_w5'] = ev(R4, 5)
log(f"R4^(5)(1/2) done ({time.time()-t0:.0f}s)")
for w in (4, 3, 2):
    vals[f'R4_w{w}'] = ev(R4, w)
vals['R4_const'] = ev(R4, 0)
vals['Psi2'] = ev(P2)
vals['Psi'] = ev(PS)
at = mp.atanh(c)
F1 = (7*c/16 - mp.mpf(5)/(8*c) + mp.mpf(3)/(16*c**3))*at + mp.mpf(9)/16 - mp.mpf(3)/(16*c**2)
F2 = (-c/16 - mp.mpf(1)/(8*c) + mp.mpf(3)/(16*c**3))*at + mp.mpf(1)/16 - mp.mpf(3)/(16*c**2)
vals['F1'] = F1; vals['F2'] = F2
z2, z3, z4 = mp.zeta(2), mp.zeta(3), mp.zeta(4)
vals['M4_hex'] = (vals['R4_w5'] + vals['R4_w4'] + vals['R4_w3'] + vals['R4_w2'] + vals['R4_const']
                  + z2*vals['Psi2'] - z3*vals['Psi'] + z4*F1 + z2**2*F2)
# weight-5 motivic part of M4 (R4^5 + zeta/log products of weight 5)
vals['M4_hex_w5part'] = vals['R4_w5']   # product parts are separate, classical by construction
log(f"B4 parts done ({time.time()-t0:.0f}s)")

# ---- Li32 atoms at {2,3}-units: single series with incremental inner sum
def li32(x, y, digits=DPS):
    x = mp.mpf(x); y = mp.mpf(y)
    r = abs(x * y)
    assert r < 1
    Nl = int((digits + 30) * mp.log(10) / -mp.log(r)) + 80
    S = mp.mpf(0)      # S_y(n) = sum_{m<n} y^m/m^2  (incremental)
    xn = x             # x^n
    ym = y             # y^m
    tot = mp.mpf(0)
    for n in range(2, Nl + 1):
        S += ym / ( (n-1) * (n-1) )   # m = n-1 term
        ym *= y
        xn *= x
        tot += xn * S / (n * n * n)
    return tot

sols = [Fr(-8), Fr(-3), Fr(-2), Fr(-1), Fr(-1,2), Fr(-1,3), Fr(-1,8), Fr(1,9), Fr(1,4),
        Fr(1,3), Fr(1,2), Fr(2,3), Fr(3,4), Fr(8,9), Fr(9,8), Fr(4,3), Fr(3,2), Fr(2), Fr(3),
        Fr(4), Fr(9), Fr(1), Fr(-1)]
xs = [s for s in sols if abs(s) <= 1 and s != 1]
ys = [s for s in sols if s != 1] + [Fr(1)]
li32_vals = {}
for x in xs:
    for y in ys:
        if abs(x * y) <= Fr(3, 4) and x != 1:
            key = f"Li32({x},{y})"
            li32_vals[key] = li32(mp.mpf(x.numerator)/x.denominator, mp.mpf(y.numerator)/y.denominator)
log(f"Li32 atoms: {len(li32_vals)} values ({time.time()-t0:.0f}s)")

# ---- reference values: Li5 basis atoms (task list), Li3 atoms, logs, zetas
def Li(s, x): return mp.polylog(s, mp.mpf(x))
li5_args = [Fr(-2), Fr(1,3), Fr(2,3), Fr(1,2), Fr(1,9), Fr(-1,2), Fr(-1,3), Fr(1,4), Fr(3,4), Fr(-1,8), Fr(8,9)]
li3_args = [Fr(1,3), Fr(-1,2), Fr(2,3), Fr(1,2), Fr(-1,3), Fr(1,4), Fr(3,4), Fr(-1,8), Fr(1,9), Fr(8,9), Fr(-2), Fr(-1)]
ref = {}
for a in li5_args: ref[f"Li5({a})"] = Li(5, a.numerator/a.denominator)
for a in li3_args: ref[f"Li3({a})"] = Li(3, a.numerator/a.denominator)
ref['Li2(1/3)'] = Li(2, mp.mpf(1)/3)
ref['Li2(2/3)'] = Li(2, mp.mpf(2)/3)
ref['Li2(1/2)'] = Li(2, mp.mpf(1)/2)
ref['Li4(1/3)'] = Li(4, mp.mpf(1)/3)
ref['Li4(2/3)'] = Li(4, mp.mpf(2)/3)
ref['Li4(1/2)'] = Li(4, mp.mpf(1)/2)
ref['Li4(1/9)'] = Li(4, mp.mpf(1)/9)
ref['Li4(3/4)'] = Li(4, mp.mpf(3)/4)
log(f"reference values done ({time.time()-t0:.0f}s)")

out = {k: mp.nstr(v, DPS) for k, v in vals.items()}
out.update({k: mp.nstr(v, DPS) for k, v in li32_vals.items()})
out.update({k: mp.nstr(v, DPS) for k, v in ref.items()})
json.dump(out, open(G2B + 'values_hex.json', 'w'), indent=1)
log(f"step_04 done ({time.time()-t0:.0f}s); saved {len(out)} values")
for k in ['R4_w5', 'M4_hex']:
    log(f"  {k} = {out[k][:60]}")
