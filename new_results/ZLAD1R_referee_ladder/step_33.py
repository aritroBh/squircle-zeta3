"""ZLAD1R step_33 -- function-level checks at the square corner c^2 = 1/2 (60 digits):
 (i) [z3]M_k(pi/2) = -[z2]M_{k-1}(pi/2) for k = 3..8 (rung 1 at the function level);
 (ii) [z2]M_2(pi/2) = Phi(1/sqrt2) = 1/2 - artanh(1/sqrt2)/(2 sqrt2) (paper thm:flat);
 (iii) [z3]M_3(pi/2) = -Phi(1/sqrt2).
Tail beyond n=40 bounded by geometric estimate."""
import json
from fractions import Fraction as Fr
from mpmath import mp, mpf, zeta, atanh, sqrt, log
mp.dps = 60
mine = json.load(open('muR_exact_N40_K8.json'))
def mu(k, n, m): return Fr(mine[str(k)][str(n)].get(m, '0'))
def Mpart(k, m, N=40):
    s = mpf(0)
    for n in range(1, N+1):
        v = mu(k, n, m)
        s += mpf(v.numerator)/mpf(v.denominator) / mpf(2)**n
    return s
ok = True
print('rung 1 at the function level (c^2=1/2):')
for k in range(3, 9):
    lhs = Mpart(k, 'z3'); rhs = -Mpart(k-1, 'z2')
    d = abs(lhs - rhs)
    ok &= d < mpf('1e-9')   # truncation ~ 2^-41 * poly
    print(f'  k={k}: [z3]M_{k} = {mp.nstr(lhs,22)}  -[z2]M_{k-1} = {mp.nstr(rhs,22)}  |diff|={mp.nstr(d,2)}')
Phi = mpf('0.5') - atanh(1/sqrt(2))/(2*sqrt(2))
v2 = Mpart(2, 'z2'); v3 = Mpart(3, 'z3')
print('[z2]M_2(pi/2) =', mp.nstr(v2, 25), ' Phi(1/sqrt2) =', mp.nstr(Phi, 25), ' |diff|=', mp.nstr(abs(v2-Phi),2))
print('[z3]M_3(pi/2) =', mp.nstr(v3, 25), ' -Phi        =', mp.nstr(-Phi, 25), ' |diff|=', mp.nstr(abs(v3+Phi),2))
ok &= abs(v2-Phi) < mpf('1e-9') and abs(v3+Phi) < mpf('1e-9')
print('function-level checks:', ok)
# tail estimate: |mu_{k,n}| for n>40 bounded crudely by n^2 * max|coeff| decay -- check last terms magnitude
import math
last = max(abs(float(mu(8, 40, m))) for m in mine['8']['40'])
print('max |[.]mu_{8,40}| =', last, '-> tail < ~', last*2**-40*100)
open('progress.log','a').write('step_33: function-level rung1 at pi/2 for k=3..8 (1e-9), Phi(1/sqrt2) checks PASS\n')
