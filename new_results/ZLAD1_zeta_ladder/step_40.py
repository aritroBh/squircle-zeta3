code = r'''"""verify_zlad1.py -- self-contained verification of the ZLAD1 zeta-ladder results.
Runs in < 5 minutes. Recomputes the exact mu_{k,n} table (k<=6, n<=80), then checks:
  (1) all ladder rungs at k<=6 for n<=80 (EXACT);
  (2) proof pillars: P2 ((d_a+d_b)E_n = u E_n) and P1 (operator formula) at sample points;
  (3) closed forms for [z5]mu_5, [z3^2]mu_6, [z2z3]mu_5, [z6]mu_6, [z2z4]mu_6, [z2^3]mu_6;
  (4) the resummation identity G_n(u) = sum_e T_e(u) Phi_e(u;zeta) numerically (40 digits, n=3).
Prints one PASS/FAIL line per group; exits nonzero on any failure.
"""
import sys, json, time
from fractions import Fraction as Fr
from math import comb, factorial

sys.path.insert(0, '.')
t0 = time.time()
N = 80
from step_02 import E_series, make_DB, poch_ratio, mu_all, const, padd, pmul, pscale, H, ZV

print('recomputing exact table (k<=6, n<=80)...', flush=True)
table = {k: {} for k in range(7)}
for nv in range(1, N+1):
    Ek = E_series(nv, 6)
    DB = make_DB(nv, 6)
    row = mu_all(nv, 6, Ek, DB)
    for k in range(7):
        table[k][nv] = row[k]
print(f'  done ({time.time()-t0:.0f}s)', flush=True)

def mu(k, n, m): return Fr(table[k][n].get(m, '0'))
fails = []
def group(name, ok):
    print(('PASS' if ok else 'FAIL'), name, flush=True)
    if not ok: fails.append(name)

# (1) rungs
ok = True
for k in range(2, 6):
    ok &= all(mu(k+1, n, 'z3') == -mu(k, n, 'z2') for n in range(1, N+1))
for k in range(4, 7):
    ok &= all(mu(k, n, 'z4') == mu(k-2, n, 'z2') + mu(k, n, 'z2^2') for n in range(1, N+1))
for k in (5, 6):
    ok &= all(mu(k, n, 'z5') == -mu(k-3, n, 'z2') - 2*mu(k-1, n, 'z2^2') for n in range(1, N+1))
ok &= all(mu(6, n, 'z6') == mu(2, n, 'z2') + 3*mu(4, n, 'z2^2') + 2*mu(6, n, 'z2^3') for n in range(1, N+1))
for k in (5, 6):
    ok &= all(mu(k, n, 'z2z3') == -2*mu(k-1, n, 'z2^2') for n in range(1, N+1))
ok &= all(mu(6, n, 'z3^2') == mu(4, n, 'z2^2') for n in range(1, N+1))
ok &= all(mu(6, n, 'z2z4') == 2*mu(4, n, 'z2^2') + 3*mu(6, n, 'z2^3') for n in range(1, N+1))
group('(1) all ladder rungs, n<=80', ok)

# (2) pillars
from step_05 import P2_check, P1_check
ok2 = P2_check()
ok2 = P1_check(2, 4) and ok2
group('(2) pillars P1, P2', ok2)

# (3) closed forms
import sympy as sp
n_ = sp.symbols('n')
one = 1/((2*n_-1)*(2*n_+1))
z22 = -n_/((2*n_-1)*(2*n_+1)*(2*n_+3))
z23 = n_*(5*n_-1)/(3*(2*n_-1)*(2*n_+1)*(2*n_+3)*(2*n_+5))
forms = [(5, 'z5', -one-2*z22), (6, 'z3^2', z22), (5, 'z2z3', -2*z22),
         (6, 'z6', one+3*z22+2*z23), (6, 'z2z4', 2*z22+3*z23), (6, 'z2^3', z23),
         (2, 'z2', one), (4, 'z2^2', z22), (4, 'z4', one+z22)]
ok3 = True
for k, m, f in forms:
    for nn in range(1, N+1):
        fv = f.subs(n_, nn)
        if Fr(int(fv.p), int(fv.q)) != mu(k, nn, m):
            ok3 = False
group('(3) closed forms', ok3)

# (4) resummation numeric check (n=3, 40 digits)
from mpmath import mp, mpf, zeta
mp.dps = 50
nv = 3; pref = -nv*(2*nv-1)
def mpfr(x): return mpf(x.numerator)/mpf(x.denominator)
def g(k, m): return mpfr(mu(k, nv, m))*pref
K = 6
G = [g(k, '1') + g(k, 'z2')*zeta(2) + g(k, 'z3')*zeta(3) + g(k, 'z4')*zeta(4)
     + g(k, 'z2^2')*zeta(2)**2 + g(k, 'z2z3')*zeta(2)*zeta(3) + g(k, 'z5')*zeta(5)
     + g(k, 'z2^3')*zeta(2)**3 + g(k, 'z2z4')*zeta(2)*zeta(4) + g(k, 'z3^2')*zeta(3)**2
     + g(k, 'z6')*zeta(6) for k in range(K+1)]
T = {e: [g(k, '1' if e == 0 else ('z2' if e == 1 else 'z2^%d'%e)) for k in range(K+1)] for e in range(4)}
u_, s_ = sp.symbols('u s')
zs = sp.symbols('z2 z3 z4 z5 z6')
p = {0: sp.Integer(2), 1: u_}
for w in range(2, 9):
    p[w] = sp.expand(u_*p[w-1] - s_*p[w-2])
logPhi = sp.expand(sum((-1)**(w+1)*zs[w-2]/w*(u_**w - p[w]) for w in range(2, 7)))
Phi_ser = sp.Integer(1); Xp = sp.Integer(1)
for j in range(1, 4):
    Xp = sp.expand(Xp*logPhi)
    Phi_ser = sp.expand(Phi_ser + Xp/factorial(j))
Phi = [sp.expand((-1)**e*factorial(e)*sp.expand(Phi_ser).coeff(s_, e)) for e in range(4)]
def numpoly(P):
    Pz = P
    for i, z in enumerate(zs):
        Pz = Pz.subs(z, zeta(i+2))
    return [mpf(sp.N(sp.expand(Pz).coeff(u_, k), 50)) if sp.expand(Pz).coeff(u_, k) != 0 else mpf(0) for k in range(K+1)]
PhiN = [numpoly(P) for P in Phi]
rec = [mpf(0)]*(K+1)
for e in range(4):
    for k1 in range(K+1):
        for k2 in range(K+1-k1):
            rec[k1+k2] += T[e][k1]*PhiN[e][k2]
maxd = max(abs(G[k]-rec[k]) for k in range(K+1))
print('  resummation max diff:', mp.nstr(maxd, 3))
group('(4) resummation G_n = sum T_e Phi_e (40 digits)', maxd < mpf('1e-40'))

print(f'total time {time.time()-t0:.0f}s')
if fails:
    print('FAILURES:', fails); sys.exit(1)
print('ALL VERIFICATIONS PASSED')
'''
open('verify_zlad1.py','w').write(code)
import subprocess, os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','verify_zlad1.py'],
                   capture_output=True, text=True, timeout=600)
print(r.stdout); print('ERR', r.stderr[-1500:])
