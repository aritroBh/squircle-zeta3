"""ZLAD1R step_24 -- end-to-end NUMERIC check: Cauchy-extract [u^k] of the actual Beta integral
G_n(u) = B(n,n)^{-1} int_0^1 y^{n-1}(1-y)^{n-1} e^{uA} (1-y+y e^{u l})^{1-2n} dy
and compare with my exact table values (zeta monomials substituted numerically)."""
import json
from mpmath import mp, mpf, mpc, pi, cos, sin, exp, log, zeta, quad
mp.dps = 60

mine = json.load(open('muR_exact_N40_K8.json'))
ZET = {f'z{w}': zeta(w) for w in range(2, 9)}
def mono_val(m):
    if m == '1': return mpf(1)
    v = mpf(1)
    i = 0
    for part in m.split('z'):
        if not part: continue
        if '^' in part:
            w, e = part.split('^'); v *= ZET['z'+w]**int(e)
        else:
            v *= ZET['z'+part]
    return v

def exact_num(k, n):
    return sum(mpf(num)/mpf(den.split('/')[-1]) if False else mpf(mpf(num.split('/')[0])/mpf(num.split('/')[1]) if '/' in num else mpf(num))*mono_val(m)
               for m, num in mine[str(k)][str(n)].items())

def Gn_integrand(y, n, u):
    if y == 0 or y == 1: return mpf(0)
    a = log(y); b = log(1-y)
    A = -(n-1)*a + n*b
    l = b - a
    f = 1 - y + y*exp(u*l)
    return y**(n-1) * (1-y)**(n-1) * exp(u*A) * f**(1-2*n)

from math import gamma as gg
def Gn(u, n):
    B = mpf(gg(n))**2/gg(2*n)
    return quad(lambda y: Gn_integrand(y, n, u), [0, 1]) / B

def coeff_k(n, k, rho=mpf('0.4')):
    f = lambda th: Gn(rho*exp(mpc(0, th)), n) * exp(mpc(0, -k*th))
    return quad(f, [0, 2*pi]) / (2*pi*rho**k)

for (k, n) in [(8, 3), (7, 5), (5, 7)]:
    mu_num = -coeff_k(n, k) / (n*(2*n-1))
    mu_ex  = exact_num(k, n)
    d = abs(mu_num - mu_ex)
    print(f'(k,n)=({k},{n}): numeric={mp.nstr(mu_num,25)}')
    print(f'           exact ={mp.nstr(mu_ex,25)}  |diff|={mp.nstr(d,3)}')
    assert d < mpf('1e-30'), 'MISMATCH'
print('NUMERIC QUADRATURE vs EXACT TABLE: ALL MATCH (>=30 digits)')
open('progress.log','a').write('step_24: Cauchy+quadrature of the Beta integral matches exact table at (8,3),(7,5),(5,7) to >=30 digits\n')
