import json, re
from fractions import Fraction
import mpmath as mp
mp.dps=120
sq2=mp.sqrt(2)
def mpq(fr): return mp.mpf(fr.numerator)/fr.denominator
def fval(pq): return mpq(pq[0])+mpq(pq[1])*sq2
def li_series(k,x):
    N=int(450/float(-mp.log10(abs(x))))+60
    s,xp=mp.mpf(0),x
    for j in range(1,N+1):
        s+=xp/mp.mpf(j)**k; xp*=x
    return s
def li32(x,y):
    N=int(450/float(-mp.log10(abs(x))))+60
    s,xn,yj,h=mp.mpf(0),x,mp.mpf(1),mp.mpf(0)
    for n in range(2,N+1):
        yj*=y; h+=yj/mp.mpf(n-1)**2; xn*=x
        s+=xn/mp.mpf(n)**3*h
    return s
args={'b':(-1,1),'-b':(1,-1),'-b2':(-3,2),'1/sq2':(0,Fraction(1,2)),'1/2':(Fraction(1,2),0),
 '1-1/sq2':(1,Fraction(-1,2)),'2-sq2':(2,-1),'2b':(-2,2),'(2-sq2)/4':(Fraction(1,2),Fraction(-1,4)),
 '(2+sq2)/4':(Fraction(1,2),Fraction(1,4)),'-b/2':(Fraction(1,2),Fraction(-1,2)),'(4-3sq2)/8':(Fraction(1,2),Fraction(-3,8))}
for nm,pq in args.items():
    x=fval(tuple(map(Fraction,pq)))
    for k in (2,5):
        d=li_series(k,x)-mp.polylog(k,x)
        if abs(d)>mp.mpf('1e-100'): print('MISMATCH li',k,nm,mp.nstr(d,3))
print('li_series all OK vs polylog')
# li32 vs nsum
for x,y in ((sq2-1,1),(1/sq2,1),(1-sq2,1),(sq2-1,-1)):
    v=li32(x,y)
    ref=mp.nsum(lambda n,m: x**n*y**m/(mp.mpf(n)**3*mp.mpf(m)**2), [2,mp.inf],[1,lambda n: n-1])
    print('li32',mp.nstr(x,4),mp.nstr(y,2),'diff',mp.nstr(v-ref,3))
