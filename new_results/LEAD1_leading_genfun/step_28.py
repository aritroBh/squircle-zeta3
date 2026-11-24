import os, json, mpmath
os.chdir('<project>/research/08_theory/LEAD1_leading_genfun')
mpmath.mp.dps = 60
from mpmath import mpf, sqrt, log, quad, pi, gamma
a = 1+sqrt(2); b = sqrt(2)-1
POLES = {'-a': -a, '0': mpf(0), '-b': -b, 'a': a}
def t_of_v(v): return (v-b)*(v+a)/(4*v)
def u_of_v(v): return (a-v)*(v+b)/(4*v)
def G(pole, x, y):
    r = POLES[pole]
    return quad(lambda v: t_of_v(v)**x * u_of_v(v)**y/(v-r), [b,1,a])
Dtop = json.load(open('Dtop.json'))
from fractions import Fraction as Fr
def S_d(d, x, y):
    tot = mpf(0)
    for name, rec in Dtop.items():
        dd, key = name.split('|'); i, j = map(int, key.split(','))
        if int(dd) != d: continue
        for pole, pr in rec.items():
            cf = Fr(pr['Dtop'][0])
            if cf: tot += mpf(cf.numerator)/cf.denominator * G(pole, x, y)
    return tot
# zigzag delta_d
import sympy as sp
X = sp.symbols('X')
ser = sp.series((1+sp.sinh(X))/sp.cosh(X), X, 0, 9).removeO()
delta = {d: mpf(str(sp.Rational(ser.coeff(X,d)))) for d in range(1,9)}
pts = [(mpf('0.4'), mpf('0.9')), (mpf('1.3'), mpf('0.6')), (mpf('0.7'), mpf('2.2')), (mpf('2.0'), mpf('2.0'))]
for d in [1,2,3]:
    for (x,y) in pts:
        s = S_d(d,x,y)
        B = gamma(x)*gamma(y)/gamma(x+y)
        print(f"d={d} (x,y)=({x},{y}): S={mpmath.nstr(s,20)}  2S/(delta_d B)={mpmath.nstr(2*s/(delta[d]*B),20)}")
