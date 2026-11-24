code = r'''"""step_03.py -- verify ALL ladder rungs exactly for n<=150 from mu_exact_N150_K6.json,
and fit closed forms for the free directions [z2]mu_{k,n}, [z2^2]mu_{k,n}, [z2^3]mu_{6,n}."""
import json
from fractions import Fraction as Fr

t = json.load(open('mu_exact_N150_K6.json'))
N = 150
def mu(k, n, m): return Fr(t[str(k)][str(n)].get(m, '0'))

checks = []
def chk(name, k, m, rhs):
    ok = all(mu(k, n, m) == rhs(n) for n in range(1, N+1))
    checks.append((name, ok))
    print(('PASS' if ok else 'FAIL'), name, flush=True)

# Rung 1: [z3]mu_{k+1} = -[z2]mu_k, k=2..5
for k in range(2, 6):
    chk(f'[z3]mu_{k+1} = -[z2]mu_{k}', k+1, 'z3', (lambda k: lambda n: -mu(k, n, 'z2'))(k))
# w=4: [z4]mu_k = [z2]mu_{k-2} + [z2^2]mu_k, k=4..6
for k in range(4, 7):
    chk(f'[z4]mu_{k} = [z2]mu_{k-2} + [z2^2]mu_{k}', k, 'z4', (lambda k: lambda n: mu(k-2, n, 'z2') + mu(k, n, 'z2^2'))(k))
# w=5: [z5]mu_k = -[z2]mu_{k-3} - 2[z2^2]mu_{k-1}, k=5,6
for k in (5, 6):
    chk(f'[z5]mu_{k} = -[z2]mu_{k-3} - 2[z2^2]mu_{k-1}', k, 'z5', (lambda k: lambda n: -mu(k-3, n, 'z2') - 2*mu(k-1, n, 'z2^2'))(k))
# w=6: [z6]mu_6 = [z2]mu_2 + 3[z2^2]mu_4 + 2[z2^3]mu_6
chk('[z6]mu_6 = [z2]mu_2 + 3[z2^2]mu_4 + 2[z2^3]mu_6', 6, 'z6', lambda n: mu(2, n, 'z2') + 3*mu(4, n, 'z2^2') + 2*mu(6, n, 'z2^3'))
# mixed: [z2z3]mu_k = -2[z2^2]mu_{k-1}, k=5,6
for k in (5, 6):
    chk(f'[z2z3]mu_{k} = -2[z2^2]mu_{k-1}', k, 'z2z3', (lambda k: lambda n: -2*mu(k-1, n, 'z2^2'))(k))
# [z3^2]mu_6 = [z2^2]mu_4
chk('[z3^2]mu_6 = [z2^2]mu_4', 6, 'z3^2', lambda n: mu(4, n, 'z2^2'))
# [z2z4]mu_6 = 2[z2^2]mu_4 + 3[z2^3]mu_6
chk('[z2z4]mu_6 = 2[z2^2]mu_4 + 3[z2^3]mu_6', 6, 'z2z4', lambda n: 2*mu(4, n, 'z2^2') + 3*mu(6, n, 'z2^3'))

print('ALL PASS:', all(c[1] for c in checks))

# ---- free directions: fit [z2]mu_{k,n} (is it rational in n? does it involve H_{n-1}?) ----
import sympy as sp
n_ = sp.symbols('n')
def fit_rational(vals, dnum, dden):
    # vals: list of Fraction for n=1..; fit (poly dnum)/(poly dden) with integer-ish coeffs
    m = dnum + dden + 1
    xs = list(range(1, m+1))
    M = []
    rhs = []
    for i, nn in enumerate(xs):
        row = [nn**j for j in range(dnum, -1, -1)] + [-(vals[nn-1])*nn**j for j in range(dden, -1, -1)]
        M.append([sp.Rational(v.p, v.q) for v in row])
        rhs.append(sp.Rational(0))
    A = sp.Matrix(M)
    ns = A.nullspace()
    if len(ns) == 1:
        v = ns[0]
        num = sum(v[j]*n_**(dnum-j) for j in range(dnum+1))
        den = sum(v[dnum+1+j]*n_**(dden-j) for j in range(dden+1))
        f = sp.simplify(num/den)
        if all(sp.Rational(f.subs(n_, nn)) == sp.Rational(vals[nn-1].p, vals[nn-1].q) for nn in range(1, N+1)):
            return f
    return None

print()
print('--- [z2]mu_{k,n} ---')
for k in range(2, 7):
    vals = [mu(k, nn, 'z2') for nn in range(1, N+1)]
    found = None
    for dd in range(0, 9):
        f = fit_rational(vals, dd+2, dd) if dd>0 else fit_rational(vals, 2, 0)
        if f is not None:
            found = f; break
    print(f'[z2]mu_{k},n =', found)
print('--- [z2^2]mu_{k,n} ---')
for k in range(4, 7):
    vals = [mu(k, nn, 'z2^2') for nn in range(1, N+1)]
    found = None
    for dd in range(0, 9):
        f = fit_rational(vals, dd+2, dd)
        if f is not None:
            found = f; break
    print(f'[z2^2]mu_{k},n =', found)
print('--- [z2^3]mu_{6,n} ---')
vals = [mu(6, nn, 'z2^3') for nn in range(1, N+1)]
found = None
for dd in range(0, 12):
    f = fit_rational(vals, dd+2, dd)
    if f is not None:
        found = f; break
print('[z2^3]mu_6,n =', found)
'''
open('step_03.py','w').write(code)
import subprocess, os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_03.py'],
                   capture_output=True, text=True, timeout=1800)
print(r.stdout)
print('ERR', r.stderr[-2000:])
