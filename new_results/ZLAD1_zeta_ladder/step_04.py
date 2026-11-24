"""step_04.py -- fit free directions [z2]mu_{k,n}, [z2^2]mu_{k,n}, [z2^3]mu_{6,n}.
Also test whether they involve harmonic numbers H_{n-1} by checking rational-in-n fits."""
import json
from fractions import Fraction as Fr
import sympy as sp

t = json.load(open('mu_exact_N150_K6.json'))
N = 150
def mu(k, n, m): return Fr(t[str(k)][str(n)].get(m, '0'))
n_ = sp.symbols('n')

def fit_rational(vals, dnum, dden):
    m = dnum + dden + 1
    M = []
    for nn in range(1, m+1):
        v = sp.Rational(vals[nn-1].numerator, vals[nn-1].denominator)
        row = [sp.Rational(nn)**j for j in range(dnum, -1, -1)]
        row += [-v*sp.Rational(nn)**j for j in range(dden, -1, -1)]
        M.append(row)
    ns = sp.Matrix(M).nullspace()
    if len(ns) == 1:
        v = ns[0]
        num = sum(v[j]*n_**(dnum-j) for j in range(dnum+1))
        den = sum(v[dnum+1+j]*n_**(dden-j) for j in range(dden+1))
        f = sp.simplify(num/den)
        if all(sp.Rational(f.subs(n_, nn)) == sp.Rational(vals[nn-1].numerator, vals[nn-1].denominator) for nn in range(1, N+1)):
            return sp.factor(f)
    return None

def try_fits(name, vals, maxd=10):
    for dd in range(0, maxd+1):
        for dn in range(dd, maxd+2):
            f = fit_rational(vals, dn, dd)
            if f is not None:
                print(name, '=', f, f'  (deg {dn}/{dd})')
                return f
    print(name, ': no rational fit up to degrees', maxd)
    return None

print('--- [z2]mu_{k,n} ---')
for k in range(2, 7):
    try_fits(f'[z2]mu_{k},n', [mu(k, nn, 'z2') for nn in range(1, N+1)])
print('--- [z2^2]mu_{k,n} ---')
for k in range(4, 7):
    try_fits(f'[z2^2]mu_{k},n', [mu(k, nn, 'z2^2') for nn in range(1, N+1)])
print('--- [z2^3]mu_{6,n} ---')
try_fits('[z2^3]mu_6,n', [mu(6, nn, 'z2^3') for nn in range(1, N+1)])

# do [z2]mu_{3,n} etc involve H_{n-1}? test: v_n = a_n + b_n H_{n-1} with a,b rational low degree
def H(m): return Fr(sum(Fr(1,i) for i in range(1, m+1))) if m >= 1 else Fr(0)
print()
print('--- H_{n-1} content test for [z2]mu_{3,n}: fit v = r1(n) + r2(n) H_{n-1} ---')
vals = [mu(3, nn, 'z2') for nn in range(1, N+1)]
# least squares style exact: unknowns coeffs of r1 (deg d1/d0), r2 (deg d1/d0): solve linear system
def fit_with_H(vals, d1, d0, useH):
    # v_n = (p(n)+q(n)H_{n-1})/(r(n)), p,q deg d1, r monic deg d0
    unk = (d1+1)*2 + d0
    rows = []
    for nn in range(1, min(unk+3, N)+1):
        v = sp.Rational(vals[nn-1].numerator, vals[nn-1].denominator)
        hv = sp.Rational(H(nn-1).numerator, H(nn-1).denominator)
        row = [sp.Rational(nn)**j for j in range(d1, -1, -1)]
        row += [hv*sp.Rational(nn)**j for j in range(d1, -1, -1)]
        row += [-v*sp.Rational(nn)**j for j in range(d0, 0, -1)]
        rows.append(row)
    ns = sp.Matrix(rows).nullspace()
    results = []
    for vec in ns:
        p = sum(vec[j]*n_**(d1-j) for j in range(d1+1))
        q = sum(vec[d1+1+j]*n_**(d1-j) for j in range(d1+1))
        r = n_**d0 + sum(vec[2*d1+2+j]*n_**(d0-1-j) for j in range(d0))
        f = sp.simplify((p + (sp.symbols('H') if useH else 0)*q)/r)
        ok = True
        for nn in range(1, N+1):
            fv = (p.subs(n_, nn) + sp.Rational(H(nn-1).numerator, H(nn-1).denominator)*q.subs(n_, nn))/r.subs(n_, nn)
            if sp.Rational(fv) != sp.Rational(vals[nn-1].numerator, vals[nn-1].denominator):
                ok = False; break
        if ok:
            results.append((p, q, r))
    return results

for k in [3, 4, 5]:
    vals = [mu(k, nn, 'z2') for nn in range(1, N+1)]
    res = None
    for d0 in range(0, 6):
        for d1 in range(0, 8):
            rr = fit_with_H(vals, d1, d0, True)
            if rr:
                res = rr[0]; break
        if res: break
    if res:
        p, q, r = res
        print(f'[z2]mu_{k},n = ({sp.expand(p)} + ({sp.expand(q)})*H_{{n-1}})/({sp.expand(r)})')
    else:
        print(f'[z2]mu_{k},n: no (p+qH)/r fit found')
