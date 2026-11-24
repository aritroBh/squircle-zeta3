"""Step 1 numerics: each form + N5 vs certified W5 (5400-bit value), and R14-R5, at DPS digits."""
import sys, json, mpmath as mp
from fractions import Fraction as F
from common import *
DPS=int(sys.argv[1]) if len(sys.argv)>1 else 360
mp.mp.dps=DPS+30
r2=mp.sqrt(2)
def fv(pq): return mp.mpf(pq[0].numerator)/pq[0].denominator + (mp.mpf(pq[1].numerator)/pq[1].denominator)*r2
cache={}
def li32(x,y):
    # sum_{n>m>0} x^n y^m/(n^3 m^2)
    s=mp.mpf(0); inner=mp.mpf(0); xn=mp.mpf(1); ym=mp.mpf(1); n=1
    eps=mp.mpf(10)**(-(DPS+20))
    while True:
        xn*=x
        term=xn*inner/mp.mpf(n)**3
        s+=term
        if n>20 and abs(term)<eps and abs(xn)<eps: break
        ym*=y; inner+=ym/mp.mpf(n)**2; n+=1
    return s
def val(k):
    if k in cache: return cache[k]
    t=k[0]
    if t=='lam': v=mp.log(2)
    elif t=='L': v=mp.log(1+r2)
    elif t=='pi2': v=mp.pi**2
    elif t=='z3': v=mp.zeta(3)
    elif t=='z5': v=mp.zeta(5)
    elif t=='Lchi': v=(mp.zeta(5,mp.mpf(1)/8)+mp.zeta(5,mp.mpf(7)/8)-mp.zeta(5,mp.mpf(3)/8)-mp.zeta(5,mp.mpf(5)/8))/mp.mpf(8)**5
    else: v=mp.polylog(k[1],fv(k[2]))
    cache[k]=v; return v
def evalform(form):
    s=mp.mpf(0)
    for key,c in terms_keys(form).items():
        p=mp.mpf(1)
        for a in key: p*=val(a)
        s+=fv(c)*p
    return s
R14=load(R14_PATH); R5=load(R5_PATH)
n5=json.load(open(N5_PATH))['correction']
N5=mp.mpf(0)
for t in n5:
    x=fv((F(t['x'][0]),F(t['x'][1]))); y=fv((F(t['y'][0]),F(t['y'][1])))
    N5+=fv((F(t['coefficient'][0]),F(t['coefficient'][1])))*li32(x,y)
W5=mp.mpf(json.load(open(W5_PATH))['value'])
v14=evalform(R14); v5=evalform(R5)
print(f'dps={DPS}  terms R14={len(R14)} R5={len(R5)}')
print('W5 (certified file)  =', mp.nstr(W5,30))
print('N5 + R14 - W5        =', mp.nstr(N5+v14-W5,5))
print('N5 + R5  - W5        =', mp.nstr(N5+v5-W5,5))
print('R14 - R5             =', mp.nstr(v14-v5,5))
