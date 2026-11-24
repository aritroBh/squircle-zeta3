import os
os.chdir('<project>/research/08_theory/LEAD1_leading_genfun')
code = r'''"""LEAD1 verification script. Reruns, from the exact kernel data, at reduced precision (dps 240):
  A. extraction of the top kernel coefficients D_d (d=1..8) and the EXACT zigzag/binomial pattern check;
  B. quadrature of the 56 master integrals and assembly of c(w), w=2..9; controls c(2),c(3),c(4),c(6);
  C. LLL identifications of c(5)..c(9) in the Gamma-world monomial basis (D=225), compared with the
     stored dps-420 relations (two-precision rule);
  D. the generating-function law  c(w) = eps_{w-2} A_{w-1}/(2 (w-2)!) [y^w] E(y)  for all w=2..9;
  E. proof-side identities: even-d and odd-d single-integral evaluations vs (d!/2) e_{d+1},
     and the F(y) partial-fraction identity.
Prints PASS/FAIL per item. Expected runtime < 15 min (mostly quadrature).
"""
import os, sys, json, time, math
os.nice(10)
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = '<project>'
THEORY = os.path.join(BASE, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')
import mpmath
mpmath.mp.dps = 240
from mpmath import mpf, sqrt, log, quad, pi, zeta, gamma
from fractions import Fraction as Fr
from math import comb
t0 = time.time()
ok_all = True
def report(name, ok, resid=None):
    global ok_all
    ok_all &= bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"   residual={mpmath.nstr(resid,6)}" if resid is not None else ""))

# ---------- A. D_d extraction + exact pattern ----------
import pickle
kern = {}
for K in range(1, 19):
    with open(os.path.join(THEORY, f'kern_{K}.pkl'), 'rb') as f:
        kern[K] = pickle.load(f)
POLE = {(0,0):"0",(1,1):"a",(-1,-1):"-a",(-1,1):"b",(1,-1):"-b"}
def ckk(k, key, r):
    P = kern.get(k, {}).get(key)
    if P is None or r not in P['real']: return (Fr(0), Fr(0))
    z = P['real'][r]; return (z[0], z[1])
Dtop = {}
for d in range(1, 9):
    keys = sorted({key for K in kern for key in kern[K] if sum(key) == d}, key=lambda kk: -kk[0])
    for key in keys:
        poles = sorted({r for K in kern for r in kern[K].get(key, {'real': {}})['real']})
        rec = {}
        for r in poles:
            seq = [ckk(K, key, r) for K in range(d+1, 19)]
            D = []; s = list(seq)
            for m in range(1, d+1):
                D.append(s[0]); s = [(x2[0]-x1[0], x2[1]-x1[1]) for x1, x2 in zip(s, s[1:])]
            rec[POLE[r]] = {"Dtop": [D[d-1][0], D[d-1][1]], "poly_ok": all(t == (Fr(0), Fr(0)) for t in s)}
        Dtop[f"{d}|{key[0]},{key[1]}"] = rec
import sympy as sp
X = sp.symbols('X')
ser = sp.series((1+sp.sinh(X))/sp.cosh(X), X, 0, 9).removeO()
delta = {d: Fr(ser.coeff(X, d).p, ser.coeff(X, d).q) for d in range(1, 9)}
def predicted(d, i, j, pole):
    de = delta[d]
    if pole == '-a': return (-1)**(j+1)*comb(d, j)*de
    if d % 2 == 1:
        if pole == 'a': return (-1)**j*comb(d, j)*de*(Fr(1) if j < (d-1)//2 else Fr(1, 2))
        if j == (d-1)//2:
            if pole == '-b': return Fr((-1)**(j+1)*comb(d, j))*de/2
            if pole == '0': return Fr((-1)**j*comb(d, j))*de/2
    else:
        if pole == '-b' and j < d//2: return (-1)**j*comb(d, j)*de
        if pole == '0' and j == d//2: return Fr((-1)**j*comb(d, j))*de/2
    return Fr(0)
bad = []
for name, rec in Dtop.items():
    d, key = name.split('|'); i, j = map(int, key.split(',')); d = int(d)
    for pole in ['-a', '0', '-b', 'a']:
        got = rec[pole]['Dtop'][0] if pole in rec else Fr(0)
        if got != predicted(d, i, j, pole): bad.append((name, pole))
poly_ok = all(pr['poly_ok'] for rec in Dtop.values() for pr in rec.values())
report("A: D_d table (d=1..8, all keys/poles) == zigzag/binomial pattern, polynomial law exact", not bad and poly_ok)

# ---------- B. quadrature + assembly + controls ----------
a = 1+sqrt(2); b = sqrt(2)-1
POLES = {'-a': -a, '0': mpf(0), '-b': -b, 'a': a}
def integrand(v, i, j, r):
    t = (v-b)*(v+a)/(4*v); u = (a-v)*(v+b)/(4*v)
    return log(t)**i * log(u)**j / (v - r)
N = {}
jobs = []
for name, rec in Dtop.items():
    d, key = name.split('|'); i, j = map(int, key.split(','))
    for pole, pr in rec.items():
        if pr['Dtop'] != [Fr(0), Fr(0)]:
            jobs.append((i, j, pole))
for (i, j, pole) in jobs:
    N[(i, j, pole)] = quad(lambda v: integrand(v, i, j, POLES[pole]), [b, 1, a])
cw = {}
for w in range(2, 10):
    d = w-1; tot = mpf(0)
    for name, rec in Dtop.items():
        dd, key = name.split('|'); i, j = map(int, key.split(','))
        if int(dd) != d: continue
        for pole, pr in rec.items():
            cf = pr['Dtop'][0]
            if cf: tot += mpf(cf.numerator)/cf.denominator * N[(i, j, pole)]
    cw[w] = tot/mpmath.factorial(d-1)
lam = log(2)
controls = {
 2: (pi**2 + 6*lam**2)/24,
 3: -(zeta(3)/8 + pi**2*lam/24 + lam**3/12),
 4: -(pi**4/320 + pi**2*lam**2/48 + lam**4/48 + zeta(3)*lam/8),
 6: (168*lam**6 + 420*pi**2*lam**4 + 378*pi**4*lam**2 + 61*pi**6 + 5040*zeta(3)*lam**3
     + 2520*pi**2*zeta(3)*lam + 3780*zeta(3)**2 + 22680*zeta(5)*lam)/mpmath.factorial(9)}
for w, tgt in controls.items():
    report(f"B: control c({w})", abs(cw[w]-tgt) < mpf(10)**(-210), cw[w]-tgt)

# ---------- C. LLL identifications ----------
from flint import fmpz_mat, fmpz
def monoms(w):
    out = []
    def rec(rem, factors, start):
        if rem == 0: out.append(list(factors)); return
        for g in range(start, rem+1):
            if g in (1, 2) or g % 2 == 1: rec(rem-g, factors+[g], g)
    rec(w, [], 1)
    return out
def mval(m):
    v = mpf(1)
    for g in m: v *= lam if g == 1 else (pi**2 if g == 2 else zeta(g))
    return v
def lll_relation(vals, D):
    n = len(vals); S = mpmath.power(10, D)
    M = []
    for i in range(n):
        row = [0]*n + [int(mpmath.nint(vals[i]*S))]; row[i] = 1; M.append(row)
    B = fmpz_mat([[fmpz(x) for x in row] for row in M]).lll(delta=mpmath.mpf('0.99'), gram='exact')
    rel = [int(B[0, i]) for i in range(n)]
    g0 = 0
    for x in rel: g0 = math.gcd(g0, abs(x))
    return [x//g0 for x in rel], sum(mpf(rel[i])*vals[i] for i in range(n))
stored = json.load(open(os.path.join(HERE, 'lll_relations_dps400.json')))
stored6 = [168, 420, 5040, 378, 2520, 22680, 61, 3780, -362880]
for w in [5, 6, 7, 8, 9]:
    basis = monoms(w)
    vals = [mval(m) for m in basis] + [cw[w]]
    rel, resid = lll_relation(vals, 225)
    tgt = stored[str(w)]['rel'] if w != 6 else stored6
    same = rel == tgt or rel == [-x for x in tgt]
    report(f"C: LLL identification c({w}) identical to dps-420 relation (D=225)", same and abs(resid) < mpf(10)**(-200), resid)

# ---------- D. generating-function law ----------
NN = 12
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
L = [mpf(0)]*(NN+1); L[1] = lam
for n in range(2, NN+1): L[n] = (zeta(n) if n % 2 == 0 else eta(n))/n
E = [mpf(0)]*(NN+1); E[0] = mpf(1)
for n in range(1, NN+1): E[n] = sum(k*L[k]*E[n-k] for k in range(1, n+1))/n
# E(y) = Gamma(1+y/2) Gamma(1/2-y/2)/sqrt(pi): series cross-check
yy = mpf('0.37')
lhsE = sum(E[n]*yy**n for n in range(NN+1))
rhsE = gamma(1+yy/2)*gamma(mpf('0.5')-yy/2)/sqrt(pi)
report("D0: E(y) series == Gamma(1+y/2)Gamma(1/2-y/2)/sqrt(pi) at y=0.37", abs(lhsE-rhsE) < mpf(10)**(-25), lhsE-rhsE)
A = [1, 1, 1, 2, 5, 16, 61, 272, 1385, 7936, 50521]
def eps(m): return [1, -1, -1, 1][m % 4]
for w in range(2, 10):
    m = w-2
    law = eps(m)*A[m+1]/(2*mpmath.factorial(m))*E[w]
    report(f"D: law c({w}) = eps A/(2 m!) e_{w}", abs(cw[w]-law) < mpf(10)**(-210), cw[w]-law)

# ---------- E. proof-side integration identities ----------
mpmath.mp.dps = 60
def evenint(d): return quad(lambda w: (-log(w))**d*((1+w)/sqrt(1+w**2)-1)/w, [0, 1])/2
def oddint(d): return quad(lambda w: (log(w))**(d+1)*(1+w)/(1+w**2)**mpf('1.5'), [0, 1])/(2*(d+1))
for d in range(1, 9):
    tgt = mpmath.factorial(d)/2*E[d+1]
    val = evenint(d) if d % 2 == 0 else oddint(d)
    report(f"E: single-integral identity d={d}", abs(val-tgt) < mpf(10)**(-45), val-tgt)
def F(y): return quad(lambda w: w**(-y)*((1+w)/sqrt(1+w**2)-1)/w, [0, 1])
def Eg(y): return gamma(1+y/2)*gamma(mpf('0.5')-y/2)/sqrt(pi)
y = mpf('0.3')
lhs = (F(y)+F(-y))/2; rhs = (Eg(y)-Eg(-y))/(2*y)
report("E: (F(y)+F(-y))/2 == (E(y)-E(-y))/(2y) at y=0.3", abs(lhs-rhs) < mpf(10)**(-40), lhs-rhs)
print(f"\nALL: {'PASS' if ok_all else 'FAIL'}   ({time.time()-t0:.0f}s, dps 240/60)")
sys.exit(0 if ok_all else 1)
'''
open('verify_lead1.py','w').write(code)
print('written', len(code))
