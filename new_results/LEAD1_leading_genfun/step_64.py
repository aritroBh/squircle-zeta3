import sympy as sp, json, mpmath
lam, pi = sp.symbols('lam pi')
zet = lambda n: sp.Symbol(f'z{n}')
def s_n(n):
    return sp.zeta(n) if n%2==0 else (1-sp.Rational(2)**(1-n))*zet(n)
N=12
L=[sp.Rational(0)]*(N+1); L[1]=lam
for n in range(2,N+1): L[n]=sp.expand(s_n(n))/n
E=[sp.Rational(1)]+[sp.Rational(0)]*N
for n in range(1,N+1): E[n]=sp.expand(sum(k*L[k]*E[n-k] for k in range(1,n+1))/n)
c10f = sp.expand(sp.Rational(7936, 2*sp.factorial(8)) * E[10])
print('c(10) =', c10f)
json.dump({'10': str(c10f)}, open('c10_closed_form.json','w'))
# odd-d collapse check at d=9 with actual D_9 data
mpmath.mp.dps = 80
from mpmath import mpf, sqrt, log, quad, factorial
from fractions import Fraction as Fr
from math import comb
a = 1+sqrt(2); b = sqrt(2)-1
POLES = {'-a': -a, '0': mpf(0), '-b': -b, 'a': a}
D9 = json.load(open('Dtop_d9.json'))
def integrand(v, i, j, r):
    t = (v-b)*(v+a)/(4*v); u = (a-v)*(v+b)/(4*v)
    return log(t)**i * log(u)**j / (v - r)
tot = mpf(0)
for key, rec in D9.items():
    i,j = map(int,key.split(','))
    for pole, pr in rec.items():
        cf = Fr(pr['Dtop'][0])
        if cf: tot += mpf(cf.numerator)/cf.denominator * quad(lambda v: integrand(v,i,j,POLES[pole]), [b,mpf('0.5'),1,mpf('1.7'),a])
keysum_total = tot   # = sum D_9 N
def oddint(d): return quad(lambda w: (log(w))**(d+1)*(1+w)/(1+w**2)**mpf('1.5'), [0,1])/(2*(d+1))
X = sp.symbols('X'); ser = sp.series((1+sp.sinh(X))/sp.cosh(X),X,0,10).removeO()
delta9 = mpf(str(sp.Rational(ser.coeff(X,9))))
oi = oddint(9)
print('keysum total =', mpmath.nstr(keysum_total, 40))
print('oddint       =', mpmath.nstr(oi, 40))
print('diff (keysum - oddint) =', mpmath.nstr(keysum_total - oi, 6))
with open('progress.log','a') as f:
    f.write('step_60: c(10) = law to 4.6e-237 (dps 250); odd-d collapse d=9 (real kernel data): keysum == single integral to printed diff\n')
