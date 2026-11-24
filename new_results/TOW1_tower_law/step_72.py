common = '''"""Common symbolic setup for the TOW1 identity proofs.
All identities are exact rational-function identities in the symbols (d,p) over Q(i,sqrt2)
(sqrt2 appears only inside explicit constants). A remainder that simplifies to 0 is a
COMPLETE PROOF for all values of the symbols (exact rational arithmetic in sympy)."""
import sympy as sp
i = sp.I
t = { (1,0): (-1-i)/2, (2,0): (-1+i)/2, (1,1): sp.Rational(3,4)*(1-i),
      (2,1): -sp.Rational(3,4)-sp.Rational(3,4)*i, (1,2): 1+sp.Rational(3,4)*i, (2,2): 1-sp.Rational(3,4)*i }
def tc1c(n):
    return (-1)**n * t[(2,n)]   # conj(tau1_n) = (-1)^n tau2_n
def rfun(dd, m, conjg=False):
    s = -2*dd if not conjg else 2*dd
    return 2*(dd-2*m) + s*i - 4*(dd-1)*(dd-2*m)/(2*dd-3)
def R2(dd, m, conjg=False):
    s = -1 if not conjg else 1
    return sp.Rational(1,2)*( (3*dd-3)/2 + s*i*rfun(dd,m,conjg)/2 - (dd-2)**2 )
'''
open('proof_common.py','w').write(common)

e1adj = '''from proof_common import *
p = sp.Symbol('p', integer=True, positive=True)
d = 2*p-1
# E1 (depth-1 identity), adjacent key (p,p-1), normalized by rho_{p,p-1}(0).
LHS = -i*(d-2)
rA0 = i*(2*d-1)/(p+1)                       # rho_{p+1,p-1}(0)/rho_{p,p-1}(0)
rA1 = -i*(d-1)*rA0
rB0 = (2*d-1)*(2*d+1)/((p+1)*p)             # rho_{p+1,p}(0)/rho_{p,p-1}(0)
rB1 = -i*d*rB0
ehat00 = t[(1,0)]*(4*p-3)/p                 # ehat_pp(0,0)/rho = h_hat_p/rho
ehat10 = (p+1)*( t[(1,0)]*rB1/(4*p-2) + t[(1,1)]*rB0/(4*p-1) )
feed_DL1 = (p+1)*( t[(1,0)]*rA1/(2*d-2) + t[(1,1)]*rA0/(2*d-1) )
feed_DL2 = p*( t[(2,0)]*ehat10/(2*d-2) + t[(2,1)]*ehat00/(2*d-1) )
ehat00c = sp.conjugate(ehat00)
ehat10c = (p+1)*( t[(2,0)]*sp.conjugate(rB1)/(4*p-2) + (-t[(2,1)])*sp.conjugate(rB0)/(4*p-1) )
g_m1c = -( ehat00c*(-t[(2,1)])/(4*p-3) + ehat10c*t[(2,0)]/(4*p-4) )
g_m0c = -( ehat00c*t[(2,0)]/(4*p-3) )
feed_fold = p*( g_m1c - i*(4*p-6)*g_m0c )
RHS = feed_DL1 + feed_DL2 + feed_fold
print('E1-adjacent remainder:', sp.cancel(sp.together(LHS - RHS)))
'''
open('prove_E1_adjacent.py','w').write(e1adj)

e2non = '''from proof_common import *
d, p = sp.symbols('d p', integer=True, positive=True)
# E2b (depth-2 value identity), non-adjacent key (p,q), q=d-p, normalized by rho_{p,q}(0).
LHS = R2(d,p) - 1/(2*d-3)
rA0 = i*(2*d-1)/(p+1)
rB0 = -i*(2*d-1)/(d-p+1)
RHS = ( (p+1)*( t[(1,2)]*rA0/(2*d-1) + t[(1,1)]*(-i*(d-1))*rA0/(2*d-2) + t[(1,0)]*(R2(d+1,p+1)*rA0)/(2*d-3) )
      + (d-p+1)*( t[(2,2)]*rB0/(2*d-1) + t[(2,1)]*(-i*(d-1))*rB0/(2*d-2) + t[(2,0)]*(R2(d+1,p)*rB0)/(2*d-3) ) )
print('E2b-non-adjacent remainder:', sp.cancel(sp.together(LHS - RHS)))
'''
open('prove_E2b_nonadjacent.py','w').write(e2non)
print('scripts written')
import subprocess
for f in ['prove_E1_adjacent.py','prove_E2b_nonadjacent.py']:
    r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python',f], capture_output=True, text=True, timeout=1200)
    print(f, '->', r.stdout.strip(), r.stderr[-500:] if r.returncode else '')
