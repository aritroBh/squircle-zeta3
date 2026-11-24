code = r'''"""step_06.py -- verify the resummed master identity for n=3, to order u^6, NUMERICALLY (40 digits):
G_n(u) = sum_{e>=0} T_e(u) * Phi_e(u;zeta),   Phi_e = (-d/ds)^e [Gamma(1+xi)Gamma(1+eta)/Gamma(1+u)]_{s=0},
xi+eta = u, xi*eta = s ; T_e(u) = [z2^e]G_n(u) from the exact table.
"""
import json
from mpmath import mp, mpf, zeta, gamma, loggamma, psi
mp.dps = 50
from fractions import Fraction as Fr

t = json.load(open('mu_exact_N150_K6.json'))
nv = 3
pref = -nv*(2*nv-1)

def g(k, mon):
    return mpf(Fr(t[str(k)][str(nv)].get(mon, '0')))*pref

K = 6
# G_3(u) coefficients g_k (numeric, from exact polynomials)
G = [g(k, '1') + g(k, 'z2')*zeta(2) + g(k, 'z3')*zeta(3) + g(k, 'z4')*zeta(4)
     + g(k, 'z2^2')*zeta(2)**2 + g(k, 'z2z3')*zeta(2)*zeta(3) + g(k, 'z5')*zeta(5)
     + g(k, 'z2^3')*zeta(2)**3 + g(k, 'z2z4')*zeta(2)*zeta(4) + g(k, 'z3^2')*zeta(3)**2
     + g(k, 'z6')*zeta(6) for k in range(K+1)]

# T_e series
T = {e: [g(k, 'z2^%d'%e if e > 1 else 'z2') for k in range(K+1)] for e in range(0, 4)}

# Phi_e: log H = sum_{w>=2} (-1)^w z(w)/w * [(u-t)^w + t^w - u^w]; H = exp(log H);
# Phi_e = (-(1/(u-2t)) d/dt)^e H |_{t=0} as u-series.
# work with truncated power series in u with t-series; use sympy for exact symbolic manipulation
import sympy as sp
u_, t_ = sp.symbols('u t')
zs = sp.symbols('z2 z3 z4 z5 z6')
def zser(expr_u):  # truncate in u at order K
    return sp.series(expr_u, u_, 0, K+1).removeO()
logH = sum((-1)**w*zs[w-2]/w*((u_-t_)**w + t_**w - u_**w) for w in range(2, 7))
Hser = sp.expand(sp.exp(zser(logH)))  # series in u,t with zeta coeffs (both small)
# Hser as poly in t with u-series coeffs: need total order: t counts like u (since (u-t)^w+t^w-u^w ~ t*u^{w-1})
Hser = sp.expand(Hser)
# truncate: keep terms with total degree (u,t) <= K
Htr = sp.Integer(0)
for mon in sp.Poly(Hser, u_, t_).terms():
    (eu, et), co = mon
    if eu + et <= K:
        Htr += co*u_**eu*t_**et
# apply (-(1/(u-2t)) d/dt)^e at t=0
def nab(F):
    d = sp.diff(F, t_)
    inv = sp.series(1/(u_-2*t_), t_, 0, K+1).removeO()
    return sp.expand(d*inv)
Phi = []
F = Htr
for e in range(0, 4):
    val = sp.expand(F).subs(t_, 0)
    Phi.append(sp.expand((-1)**e*val))
    F = nab(F)
# numeric Phi_e as u-series
def numpoly(P):
    Pz = P
    for i, z in enumerate(zs):
        Pz = Pz.subs(z, zeta(i+2))
    return [mpf(sp.N(Pz.coeff(u_, k), 50)) for k in range(K+1)]
PhiN = [numpoly(P) for P in Phi]
# reconstruct
rec = [mpf(0)]*(K+1)
for e in range(4):
    for k1 in range(K+1):
        for k2 in range(K+1-k1):
            rec[k1+k2] += T[e][k1]*PhiN[e][k2]
print('k : G_3(u) exact        reconstruction     diff')
for k in range(K+1):
    print(k, mp.nstr(G[k], 12), mp.nstr(rec[k], 12), mp.nstr(abs(G[k]-rec[k]), 3))
'''
open('step_06.py','w').write(code)
import subprocess, os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_06.py'],
                   capture_output=True, text=True, timeout=1200)
print(r.stdout); print('ERR', r.stderr[-1500:])
