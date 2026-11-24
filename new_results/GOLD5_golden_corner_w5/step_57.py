step02 = r'''"""step_02.py -- GOLD5 Task 1: M4 and its weight-5 part W5 = R4^(5) at the golden corners
c0 = 1/phi and c0 = sqrt5-2, from B4 (G2's term lists), to ~1450 digits.

Method: planner_check_b4.py's independent mpmath HPL series evaluator (Remiddi-Vermaseren
recursion), series order N from |c0|. H_w(-c) evaluated directly (the H_e-sign bug from G2b
stays fixed: kind 'He' adds (-c)^e H_w(-c) with H_w at -c evaluated directly).
Cross-checks:
 (a) exact flat series: flat_exact_N640_K4.json summed at c0^2 (converges like c0^{2n};
     tail ~ c0^{1282} ~ 1e-268 (gol1) resp. 1e-803 (gol3));
 (b) corner quadrature (G1 corner_series.corner_coeffs) at dps ~100.
Saves values_golden.json (1450-digit strings).
"""
import sys, os, json, re, time, mpmath as mp
from fractions import Fraction as Fr
BASE = "<project>/research/08_theory"
GOLD = BASE + "/GOLD5_golden_corner_w5/"
G2 = BASE + "/G2_M4alpha/"
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
sys.path.insert(0, G1)
DPS = 1450
mp.mp.dps = DPS + 40
t0 = time.time()

def log(msg):
    print(msg, flush=True)
    with open(GOLD + "progress.log", "a") as f:
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
log(f"step_02: dps {DPS}; terms loaded ({time.time()-t0:.0f}s)")

# ---- exact flat series at c^2
flat = json.load(open(G2 + 'flat_exact_N640_K4.json'))
KEY = {'(0, 0, 0)': '1', '(1, 0, 0)': 'z2', '(0, 1, 0)': 'z3', '(0, 0, 1)': 'z4', '(2, 0, 0)': 'z2^2'}
mu = []
for nv in range(1, 641):
    acc = {}
    for k in '01234':
        for kk, v in flat[k][str(nv)].items():
            acc[KEY[kk]] = Fr(acc.get(KEY[kk], 0)) + Fr(v)
    mu.append(acc)
def flat_val(c):
    z2 = mp.pi**2/6; z3 = mp.zeta(3); z4 = mp.pi**4/90
    vals = {'1': mp.mpf(1), 'z2': z2, 'z3': z3, 'z4': z4, 'z2^2': z2**2}
    c2 = c*c; cp = c2; tot = mp.mpf(0)
    for nv in range(640):
        for kk, v in mu[nv].items():
            tot += mp.mpf(v.numerator)/v.denominator*vals[kk]*cp
        cp *= c2
    # tail bound: |c2|^641 times a generous factor
    return tot, abs(c2)**641 * mp.mpf(100)

import corner_series as cs

def evaluate(tag, c):
    N = int((DPS + 30) / -mp.log10(abs(c))) + 60
    log(f"{tag}: c = {mp.nstr(c, 20)}, series order N = {N} ({time.time()-t0:.0f}s)")
    cache = {}
    def H(w, x):
        if (w, x) not in cache: cache[(w, x)] = hpl(w, x, N)
        return cache[(w, x)]
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
                f += (-c) ** e * H(w, -c)
            tot += qf * f
        return tot
    vals = {}
    vals['R4_w5'] = ev(R4, 5)
    log(f"  R4^(5) done ({time.time()-t0:.0f}s)")
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
    vals['M4'] = (vals['R4_w5'] + vals['R4_w4'] + vals['R4_w3'] + vals['R4_w2'] + vals['R4_const']
                  + z2*vals['Psi2'] - z3*vals['Psi'] + z4*F1 + z2**2*F2)
    vals['W5'] = vals['R4_w5']
    log(f"  B4 assembled; M4 = {mp.nstr(vals['M4'], 40)} ({time.time()-t0:.0f}s)")
    # (a) flat-series cross-check
    fv, tail = flat_val(c)
    d_flat = abs(vals['M4'] - fv)
    log(f"  (a) |B4 - flat_N640| = {mp.nstr(d_flat, 3)} (flat tail bound {mp.nstr(tail, 3)})")
    vals['flat_diff'] = d_flat; vals['flat_tail'] = tail
    return vals

def quadrature_check(tag, c, vals):
    d = mp.mp.dps; mp.mp.dps = 110
    a = 2*mp.acos(c)
    qv = cs.corner_coeffs(a, 4)[4]
    mp.mp.dps = d
    dq = abs(vals['M4'] - qv)
    log(f"  (b) |B4 - quadrature(dps 110)| = {mp.nstr(dq, 3)} (a = 2*acos(c0) = {mp.nstr(a, 12)})")
    vals['quad_diff'] = dq
    vals['quad_val'] = qv

out = {}
r5 = mp.sqrt(5)
c1 = (r5 - 1)/2
c3 = r5 - 2
for tag, c in [('gol1', c1), ('gol3', c3)]:
    vals = evaluate(tag, c)
    quadrature_check(tag, c, vals)
    out[tag] = {k: mp.nstr(v, DPS) if isinstance(v, mp.mpf) else str(v) for k, v in vals.items()}
json.dump(out, open(GOLD + 'values_golden.json', 'w'), indent=1)
log(f"step_02 done ({time.time()-t0:.0f}s); saved values_golden.json")
'''
open('step_02.py', 'w').write(step02)
import subprocess
r = subprocess.run(['nice', '-n', '10', '<project>/.venv/bin/python', 'step_02.py'], capture_output=True, text=True, cwd='.', timeout=2400)
print(r.stdout)
print("STDERR:", r.stderr[-1500:])
