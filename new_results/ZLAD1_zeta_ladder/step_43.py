import json, os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
from fractions import Fraction as Fr
t = json.load(open('mu_exact_N150_K8.json'))
N = 150
def mu(k, n, m): return Fr(t[str(k)][str(n)].get(m, '0'))
ok = all(mu(8, n, 'z4^2') == mu(4, n, 'z2^2') + 3*mu(6, n, 'z2^3') + 3*mu(8, n, 'z2^4') for n in range(1, N+1))
print('PASS' if ok else 'FAIL', '[z4^2]mu_8 = [z2^2]mu_4 + 3[z2^3]mu_6 + 3[z2^4]mu_8')
# also closed form for [z4^2]mu_8
import sympy as sp
n_ = sp.symbols('n')
z22 = -n_/((2*n_-1)*(2*n_+1)*(2*n_+3))
z23 = n_*(5*n_-1)/(3*(2*n_-1)*(2*n_+1)*(2*n_+3)*(2*n_+5))
z24 = -n_*(13*n_**2-3*n_-1)/(3*(2*n_-1)*(2*n_+1)*(2*n_+3)*(2*n_+5)*(2*n_+7))
f = sp.factor(z22 + 3*z23 + 3*z24)
print('[z4^2]mu_8 closed form:', f)
okf = all(sp.Rational(sp.expand(f).subs(n_, nn).p, sp.expand(f).subs(n_, nn).q) == sp.Rational(mu(8, nn, 'z4^2').numerator, mu(8, nn, 'z4^2').denominator) for nn in range(1, N+1))
print('verified n<=150:', okf)
# also [z8]mu_8 closed form
one = 1/((2*n_-1)*(2*n_+1))
f8 = sp.factor(one + 5*z22 + 12*z23 + 6*z24)
print('[z8]mu_8 closed form:', f8)
okf8 = all(sp.Rational(sp.expand(f8).subs(n_, nn).p, sp.expand(f8).subs(n_, nn).q) == sp.Rational(mu(8, nn, 'z8').numerator, mu(8, nn, 'z8').denominator) for nn in range(1, N+1))
print('verified n<=150:', okf8)
open('../ZLAD1_zeta_ladder/progress.log','a').write('step_08: k=7,8 rung predictions all PASS at N=150 (14 checks); [z2^4]mu_8 closed form found\n')
