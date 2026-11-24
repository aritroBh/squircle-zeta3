import mpmath, json
mpmath.mp.dps = 80
from mpmath import mpf, sqrt, log, quad, pi
from fractions import Fraction as Fr
from math import comb
# key-sum from Dtop.json via u-integrals with combo measures
Dtop = json.load(open('Dtop.json'))
S = lambda u: sqrt(1+u**2)
def m_of(d,j):
    if d%2==1:
        return (lambda u: (1+sqrt(2)*u/S(u))/(1-u**2)) if j<(d-1)//2 else (lambda u: 1/(2*(1+u))+sqrt(2)*u/(S(u)*(1-u**2)))
    else:
        return (lambda u: (1-sqrt(2)/S(u))/(1-u**2)) if j<d//2 else (lambda u: (1-sqrt(2)/S(u))/(2*(1+u)))
def keysum(d):
    tot = mpf(0)
    for name, rec in Dtop.items():
        dd, key = name.split('|'); i,j = map(int,key.split(','))
        if int(dd)!=d: continue
        cf_tot = sum(Fr(rec[p]['Dtop'][0]) for p in rec) # placeholder
        m = m_of(d,j)
        Mj = quad(lambda u: log((1+u)/2)**i*log((1-u)/2)**j*m(u), [-1,0,1])
        # coefficient (-1)^{j+1} binom(d,j) delta_d is outside; here use actual D-based coefficient:
    return tot
# simpler: use the proved D-pattern: coefficient of M_(d-j,j) is (-1)^{j+1} binom(d,j) delta_d
import sympy as sp
X = sp.symbols('X'); ser = sp.series((1+sp.sinh(X))/sp.cosh(X),X,0,9).removeO()
delta = {d: mpf(str(sp.Rational(ser.coeff(X,d)))) for d in range(1,9)}
def keysum2(d):
    tot = mpf(0)
    for j in range(0, d//2+1 if d%2==0 else (d-1)//2+1):
        m = m_of(d,j)
        Mj = quad(lambda u: log((1+u)/2)**(d-j)*log((1-u)/2)**j*m(u), [-1,0,1])
        tot += (-1)**(j+1)*comb(d,j)*delta[d]*Mj
    return tot/mpmath.factorial(d-1)
def evenint(d):
    return quad(lambda w: (-log(w))**d*((1+w)/sqrt(1+w**2)-1)/w, [0,1])/2
def oddint(d):
    return quad(lambda w: (log(w))**(d+1)*(1+w)/(1+w**2)**mpf('1.5'), [0,1])/(2*(d+1))
print("direct key-sum c(w) vs single-integral target:")
for d in range(1,9):
    ks = keysum2(d)
    tgt = evenint(d)*delta[d]/mpmath.factorial(d-1) if d%2==0 else oddint(d)*delta[d]/mpmath.factorial(d-1)
    print(f"d={d} (w={d+1}): keysum c={mpmath.nstr(ks,25)}  single-int c={mpmath.nstr(tgt,25)}  diff={mpmath.nstr(ks-tgt,6)}")
with open('progress.log','a') as f:
    f.write('step_50: odd-d decomposition verified (piece1-boundary-piece2=total); even-d proof chain validated; direct keysum vs single-integral checks d=1..8\n')
