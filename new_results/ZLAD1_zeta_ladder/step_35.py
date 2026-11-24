import os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
code = r'''"""step_06.py -- verify the resummed master identity for n=3, to order u^6, NUMERICALLY (50 digits):
G_n(u) = sum_{e>=0} T_e(u) * Phi_e(u;zeta),
Phi_e = (-1)^e e! [s^e] exp( sum_{w>=2} (-1)^{w+1} zeta(w) beta_w(u,s)/w )  (= (-d_s)^e Gamma(1+xi)Gamma(1+eta)/Gamma(1+u))."""
import json
from mpmath import mp, mpf, zeta
mp.dps = 50
from fractions import Fraction as Fr
from math import factorial, comb
import sympy as sp

t = json.load(open('mu_exact_N150_K6.json'))
nv = 3
pref = -nv*(2*nv-1)
def mpfr(x): return mpf(x.numerator)/mpf(x.denominator)
def g(k, mon): return mpfr(Fr(t[str(k)][str(nv)].get(mon, '0')))*pref
K = 6
G = [g(k, '1') + g(k, 'z2')*zeta(2) + g(k, 'z3')*zeta(3) + g(k, 'z4')*zeta(4)
     + g(k, 'z2^2')*zeta(2)**2 + g(k, 'z2z3')*zeta(2)*zeta(3) + g(k, 'z5')*zeta(5)
     + g(k, 'z2^3')*zeta(2)**3 + g(k, 'z2z4')*zeta(2)*zeta(4) + g(k, 'z3^2')*zeta(3)**2
     + g(k, 'z6')*zeta(6) for k in range(K+1)]
T = {e: [g(k, 'z2^%d'%e if e > 1 else 'z2') for k in range(K+1)] for e in range(0, 4)}

u_, s_ = sp.symbols('u s')
zs = sp.symbols('z2 z3 z4 z5 z6')
# Lucas p_w and beta_w
p = {0: sp.Integer(2), 1: u_}
for w in range(2, 9):
    p[w] = sp.expand(u_*p[w-1] - s_*p[w-2])
logPhi = sp.expand(sum((-1)**(w+1)*zs[w-2]/w*(u_**w - p[w]) for w in range(2, 7)))
Phi_ser = sp.Integer(1)
Xp = sp.Integer(1)
for j in range(1, 4):
    Xp = sp.expand(Xp*logPhi)
    Phi_ser = sp.expand(Phi_ser + Xp/factorial(j))
Phi = []
for e in range(0, 4):
    ce = sp.expand(Phi_ser).coeff(s_, e)
    Phi.append(sp.expand((-1)**e*factorial(e)*ce))
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
print('k :   G_3(u) direct            reconstruction          |diff|')
maxd = mpf(0)
for k in range(K+1):
    d = abs(G[k]-rec[k]); maxd = max(maxd, d)
    print(k, mp.nstr(G[k], 15), mp.nstr(rec[k], 15), mp.nstr(d, 3))
print('MAX DIFF:', mp.nstr(maxd, 3))
# also print the symbolic Phi_1, Phi_2 for the record
print('Phi_1 =', Phi[1])
print('Phi_2 =', Phi[2])
'''
open('step_06.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_06.py'],
                   capture_output=True, text=True, timeout=1200)
print(r.stdout); print('ERR', r.stderr[-1500:])
