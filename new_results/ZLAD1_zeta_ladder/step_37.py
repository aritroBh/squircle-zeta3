code = r'''"""step_07.py -- consequences at the square corner (c^2 = 1/2), 60 digits.
Evaluate sum_n [mon]mu_{k,n} / 2^n for various monomials and k; check constancy in k;
validate the full M_k(pi/2) against known squircle values (M4 = 0.0060279931126091225045161082508776...).
"""
import json
from mpmath import mp, mpf, zeta, log, sqrt, nstr
mp.dps = 60
from fractions import Fraction as Fr

t = json.load(open('mu_exact_N150_K6.json'))
def mpfr(x): return mpf(x.numerator)/mpf(x.denominator)
def mu(k, n, m): return mpfr(Fr(t[str(k)][str(n)].get(m, '0')))

def Mk_mon(k, m, tail=True):
    # sum_{n>=1} [mon]mu_{k,n} 2^{-n}, with geometric-tail estimate from last terms
    s = mpf(0)
    last = []
    for n in range(1, 151):
        v = mu(k, n, m)*mpf(2)**(-n)
        s += v
        last.append(v)
    return s

def Mk_full(k):
    # full numeric M_k(pi/2) from all monomials
    s = mpf(0)
    for n in range(1, 151):
        row = t[str(k)][str(n)]
        v = mpfr(Fr(row.get('1', '0')))
        for m, zz in [('z2', zeta(2)), ('z3', zeta(3)), ('z4', zeta(4)), ('z2^2', zeta(2)**2),
                      ('z2z3', zeta(2)*zeta(3)), ('z5', zeta(5)), ('z2^3', zeta(2)**3),
                      ('z2z4', zeta(2)*zeta(4)), ('z3^2', zeta(3)**2), ('z6', zeta(6))]:
            v += mpfr(Fr(row.get(m, '0')))*zz
        s += v*mpf(2)**(-n)
    return s

print('--- full M_k(pi/2) vs known ---')
M4_known = mpf('0.0060279931126091225045161082508776')
print('M4 flat series:', nstr(Mk_full(4), 30))
print('M4 known      :', nstr(M4_known, 30), ' diff:', nstr(abs(Mk_full(4)-M4_known), 3))
print('M2 flat series:', nstr(Mk_full(2), 30))
print('M3 flat series:', nstr(Mk_full(3), 30))

print()
print('--- [z2] M_k^{flat}(pi/2) = sum_n [z2]mu_{k,n}/2^n ---')
for k in range(2, 7):
    print(k, nstr(Mk_mon(k, 'z2'), 40))
print()
print('--- [z3] M_k^{flat}(pi/2) ---')
for k in range(3, 7):
    print(k, nstr(Mk_mon(k, 'z3'), 40))
print()
print('--- [z2^2] M_k^{flat}(pi/2) ---')
for k in range(4, 7):
    print(k, nstr(Mk_mon(k, 'z2^2'), 40))
print()
print('--- [z3^2] M_6, [z2^3] M_6, [z4] M_4 ---')
print('[z3^2]M6:', nstr(Mk_mon(6, 'z3^2'), 40))
print('[z2^3]M6:', nstr(Mk_mon(6, 'z2^3'), 40))
print('[z4]M4  :', nstr(Mk_mon(4, 'z4'), 40))
print('[z5]M5  :', nstr(Mk_mon(5, 'z5'), 40))
print('[z6]M6  :', nstr(Mk_mon(6, 'z6'), 40))
# closed forms for reference
c = 1/sqrt(2)
Phi = mpf(1)/2 - (1-c**2)/(2*c)*(mpf(1)/2*log((1+c)/(1-c)))
print()
print('Phi(1/sqrt2) = 1/2 - atanh(1/sqrt2)/(2sqrt2) =', nstr(Phi, 40))
print('T wt2 paper: (pi^2+6 log2^2)/24; z2-coeff = 1/4 =', mpf(1)/4)
'''
open('step_07.py','w').write(code)
import subprocess, os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_07.py'],
                   capture_output=True, text=True, timeout=1800)
print(r.stdout); print('ERR', r.stderr[-1500:])
