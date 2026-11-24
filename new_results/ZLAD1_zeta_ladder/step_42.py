import json, os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
from fractions import Fraction as Fr
import sympy as sp
t = json.load(open('mu_exact_N150_K8.json'))
N = 150
def mu(k, n, m): return Fr(t[str(k)][str(n)].get(m, '0'))
checks = []
def chk(name, k, m, rhs):
    ok = all(mu(k, n, m) == rhs(n) for n in range(1, N+1))
    checks.append((name, ok))
    print(('PASS' if ok else 'FAIL'), name)

# single-zeta w=7,8
for k in (7, 8):
    chk(f'[z7]mu_{k} = -[z2]mu_{k-5} - 4[z2^2]mu_{k-3} - 6[z2^3]mu_{k-1}', k, 'z7',
        (lambda k: lambda n: -mu(k-5, n, 'z2') - 4*mu(k-3, n, 'z2^2') - 6*mu(k-1, n, 'z2^3'))(k))
chk('[z8]mu_8 = [z2]mu_2 + 5[z2^2]mu_4 + 12[z2^3]mu_6 + 6[z2^4]mu_8', 8, 'z8',
    lambda n: mu(2, n, 'z2') + 5*mu(4, n, 'z2^2') + 12*mu(6, n, 'z2^3') + 6*mu(8, n, 'z2^4'))
# mixed weight 7
for k in (7, 8):
    chk(f'[z2z5]mu_{k} = -2[z2^2]mu_{k-3} - 6[z2^3]mu_{k-1}', k, 'z2z5',
        (lambda k: lambda n: -2*mu(k-3, n, 'z2^2') - 6*mu(k-1, n, 'z2^3'))(k))
    chk(f'[z3z4]mu_{k} = -2[z2^2]mu_{k-3} - 3[z2^3]mu_{k-1}', k, 'z3z4',
        (lambda k: lambda n: -2*mu(k-3, n, 'z2^2') - 3*mu(k-1, n, 'z2^3'))(k))
    chk(f'[z2^2z3]mu_{k} = -3[z2^3]mu_{k-1}', k, 'z2^2z3',
        (lambda k: lambda n: -3*mu(k-1, n, 'z2^3'))(k))
# mixed weight 8
chk('[z2z3^2]mu_8 = 3[z2^3]mu_6', 8, 'z2z3^2', lambda n: 3*mu(6, n, 'z2^3'))
chk('[z2^2z4]mu_8 = 3[z2^3]mu_6 + 6[z2^4]mu_8', 8, 'z2^2z4', lambda n: 3*mu(6, n, 'z2^3') + 6*mu(8, n, 'z2^4'))
chk('[z3z5]mu_8 = 2[z2^2]mu_4 + 6[z2^3]mu_6', 8, 'z3z5', lambda n: 2*mu(4, n, 'z2^2') + 6*mu(6, n, 'z2^3'))
chk('[z4^2]mu_8 = [z2^2]mu_4 - 3[z2^3]mu_6 + 3[z2^4]mu_8', 8, 'z4^2', lambda n: mu(4, n, 'z2^2') - 3*mu(6, n, 'z2^3') + 3*mu(8, n, 'z2^4'))
chk('[z2z6]mu_8 = 2[z2^2]mu_4 + 9[z2^3]mu_6 + 8[z2^4]mu_8', 8, 'z2z6', lambda n: 2*mu(4, n, 'z2^2') + 9*mu(6, n, 'z2^3') + 8*mu(8, n, 'z2^4'))
print('ALL k=7,8 RUNG PREDICTIONS PASS:', all(c[1] for c in checks))

# fit [z2^4]mu_{8,n} closed form (expect deg 3 numerator / product odd denoms)
n_ = sp.symbols('n')
z24 = n_*(5*n_-1)/(3*(2*n_-1)*(2*n_+1)*(2*n_+3)*(2*n_+5))
vals = [mu(8, nn, 'z2^4') for nn in range(1, N+1)]
# try numerator poly deg 3 over (2n-1)(2n+1)(2n+3)(2n+5)(2n+7)
den = sp.prod(2*n_+2*i-1 for i in range(5))
numf = None
import itertools
sol = sp.solve([sp.Eq(sum(sp.symbols('a0:4')[j]*nn**j for j in range(4)), sp.Rational(vals[nn-1].numerator, vals[nn-1].denominator)*den.subs(n_, nn)) for nn in range(1, 6)], sp.symbols('a0:4'))
print('[z2^4]mu_8 numerator solve:', sol)
if sol:
    num = sum(sol[sp.symbols('a0:4')[j]]*n_**j for j in range(4))
    f = sp.factor(num/den)
    okf = all(sp.Rational(sp.expand(f).subs(n_, nn).p, sp.expand(f).subs(n_, nn).q) == sp.Rational(vals[nn-1].numerator, vals[nn-1].denominator) for nn in range(1, N+1))
    print('[z2^4]mu_8 =', f, ' verified n<=150:', okf)
