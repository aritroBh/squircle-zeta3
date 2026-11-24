import json, os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
from fractions import Fraction as Fr
import sympy as sp
t = json.load(open('mu_exact_N300_K6.json'))
N = 300
def mu(k, n, m): return Fr(t[str(k)][str(n)].get(m, '0'))
checks = []
def chk(name, k, m, rhs):
    ok = all(mu(k, n, m) == rhs(n) for n in range(1, N+1))
    checks.append((name, ok))
for k in range(2, 6):
    chk(f'rung1 k={k}', k+1, 'z3', (lambda k: lambda n: -mu(k, n, 'z2'))(k))
for k in range(4, 7):
    chk(f'z4 k={k}', k, 'z4', (lambda k: lambda n: mu(k-2, n, 'z2') + mu(k, n, 'z2^2'))(k))
for k in (5, 6):
    chk(f'z5 k={k}', k, 'z5', (lambda k: lambda n: -mu(k-3, n, 'z2') - 2*mu(k-1, n, 'z2^2'))(k))
chk('z6 k=6', 6, 'z6', lambda n: mu(2, n, 'z2') + 3*mu(4, n, 'z2^2') + 2*mu(6, n, 'z2^3'))
for k in (5, 6):
    chk(f'z2z3 k={k}', k, 'z2z3', (lambda k: lambda n: -2*mu(k-1, n, 'z2^2'))(k))
chk('z3^2 k=6', 6, 'z3^2', lambda n: mu(4, n, 'z2^2'))
chk('z2z4 k=6', 6, 'z2z4', lambda n: 2*mu(4, n, 'z2^2') + 3*mu(6, n, 'z2^3'))
print('ALL RUNGS PASS at N=300:', all(c[1] for c in checks), f'({len(checks)} checks)')

n_ = sp.symbols('n')
one = 1/((2*n_-1)*(2*n_+1))
z22 = -n_/((2*n_-1)*(2*n_+1)*(2*n_+3))
z23 = n_*(5*n_-1)/(3*(2*n_-1)*(2*n_+1)*(2*n_+3)*(2*n_+5))
forms = {
 '[z5]mu_5': -one - 2*z22,
 '[z3^2]mu_6': z22,
 '[z2z3]mu_5': -2*z22,
 '[z6]mu_6': one + 3*z22 + 2*z23,
 '[z2z4]mu_6': 2*z22 + 3*z23,
}
for name, f in forms.items():
    print(name, '=', sp.factor(f))
# verify against table at all n<=300
def chkform(k, m, f):
    return all(sp.Rational(sp.expand(f).subs(n_, nn).p, sp.expand(f).subs(n_, nn).q) == sp.Rational(mu(k, nn, m).numerator, mu(k, nn, m).denominator) for nn in range(1, N+1))
print('closed forms vs table (N=300):')
print(' [z5]mu_5 = -3/((2n+3)(4n^2-1)):', chkform(5, 'z5', forms['[z5]mu_5']))
print(' [z3^2]mu_6:', chkform(6, 'z3^2', forms['[z3^2]mu_6']))
print(' [z2z3]mu_5:', chkform(5, 'z2z3', forms['[z2z3]mu_5']))
print(' [z6]mu_6:', chkform(6, 'z6', forms['[z6]mu_6']))
print(' [z2z4]mu_6:', chkform(6, 'z2z4', forms['[z2z4]mu_6']))
print(' [z2^3]mu_6 = n(5n-1)/(3(2n-1)(2n+1)(2n+3)(2n+5)):', chkform(6, 'z2^3', z23))
