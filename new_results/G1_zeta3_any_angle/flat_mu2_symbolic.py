"""flat_mu2_symbolic.py -- closed formula for mu_{1,n} and mu_{2,n} (symbolic n) in terms of
   h1 = H_{n-1}, h2 = H_{2n-1}, g1 = H^(2)_{n-1}, g2 = H^(2)_{2n-1}   (and zeta(2)),
using psi(n+j) = psi(n) + sum_{t<j} 1/(n+t), psi(2n+j) likewise, and the moment formulas
  D_{1,0} = F_p, D_{0,1} = F_q, D_{2,0} = F_p^2 + F_pp, D_{1,1} = F_p F_q + F_pq, D_{0,2} = F_q^2 + F_qq,
  F_p = psi(p)-psi(p+q), F_q = psi(q)-psi(p+q), F_pp = psi'(p)-psi'(p+q), F_pq = -psi'(p+q),
  F_qq = psi'(q)-psi'(p+q),  at p = n+j, q = n.
Checks the result against flat_exact (exact rationals) for n = 1..30.
"""
import json
from fractions import Fraction as Fr
import sympy as sp

u, y, a, b, n = sp.symbols('u y a b n')
h1, h2, g1, g2, z2 = sp.symbols('h1 h2 g1 g2 z2')   # H_{n-1}, H_{2n-1}, H2_{n-1}, H2_{2n-1}
KM = 2
A = -(n - 1)*a + n*b
l = b - a
logE = u*A + (1 - 2*n)*sp.log(1 + y*(sp.exp(u*l) - 1))
E_ser = sp.series(sp.exp(sp.series(logE, u, 0, KM + 1).removeO()), u, 0, KM + 1).removeO()


def psi0(arg):          # psi(arg)+gamma for arg = n + j or 2n + j (j >= 0 integer)
    arg = sp.expand(arg)
    if arg.coeff(n) == 1:
        j = int(arg - n)
        return h1 + sum(sp.Integer(1)/(n + t) for t in range(0, j)) if j >= 0 else None
    j = int(arg - 2*n)
    return h2 + sum(sp.Integer(1)/(2*n + t) for t in range(0, j))


def psi1(arg):          # psi'(arg) = z2 - H2_{arg-1}
    arg = sp.expand(arg)
    if arg.coeff(n) == 1:
        j = int(arg - n)
        return z2 - g1 - sum(sp.Integer(1)/(n + t)**2 for t in range(0, j))
    j = int(arg - 2*n)
    return z2 - g2 - sum(sp.Integer(1)/(2*n + t)**2 for t in range(0, j))


def D(i, l_, p, q):
    Fp = psi0(p) - psi0(p + q); Fq = psi0(q) - psi0(p + q)
    Fpp = psi1(p) - psi1(p + q); Fpq = -psi1(p + q); Fqq = psi1(q) - psi1(p + q)
    return {(0, 0): 1, (1, 0): Fp, (0, 1): Fq, (2, 0): Fp**2 + Fpp, (1, 1): Fp*Fq + Fpq,
            (0, 2): Fq**2 + Fqq}[(i, l_)]


def poch(x, j):
    return sp.Mul(*[x + t for t in range(j)]) if j > 0 else sp.Integer(1)


res = {}
for k in (1, 2):
    P = sp.Poly(sp.expand(E_ser.coeff(u, k)), y, a, b)
    tot = 0
    for (j, i, ll), co in P.terms():
        tot += co*poch(n, j)/poch(2*n, j)*D(i, ll, n + j, n)
    mu = sp.expand(-tot/(n*(2*n - 1)))
    mu = sp.collect(mu, [h1, h2, g1, g2, z2], evaluate=False)
    res[k] = {str(key): sp.factor(sp.simplify(val)) for key, val in mu.items()}
    print(f'mu_{{{k},n}} =')
    for key, val in res[k].items():
        print(f'    [{key}] : {val}')

# check against exact table (pure Fraction arithmetic, all n <= 320)
tab = json.load(open('flat_exact_N320_K3.json'))
from functools import lru_cache
@lru_cache(None)
def Hs(r, m):
    return Fr(0) if m == 0 else Hs(r, m - 1) + Fr(1, m**r)
fns = {k: {key: (lambda *v, expr=expr: expr.subs(dict(zip([n, h1, h2, g1, g2], [sp.Rational(x.numerator, x.denominator) for x in v])))) for key, expr in res[k].items()} for k in (1, 2)}
ok = True
for k in (1, 2):
    for nv in range(1, 321):
        vals = (Fr(nv), Hs(1, nv - 1), Hs(1, 2*nv - 1), Hs(2, nv - 1), Hs(2, 2*nv - 1))
        env = dict(zip(['n', 'h1', 'h2', 'g1', 'g2'], vals))
        r = Fr(0); zc = Fr(0)
        for key, f in fns[k].items():
            vv = sp.Rational(f(*vals)); v = Fr(int(vv.p), int(vv.q))
            if key == 'z2':
                zc += v
            elif key == '1':
                r += v
            else:
                r += v*Fr(eval(key, {}, env))
        exact = tab[str(k)][str(nv)]
        if r != Fr(exact.get('(0, 0, 0)', '0')) or zc != Fr(exact.get('(1, 0, 0)', '0')):
            ok = False
            print('MISMATCH', k, nv)
            break
print('symbolic harmonic-number formula agrees with exact table for k=1,2 and all n<=320:', ok)
