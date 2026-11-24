import json, re
from fractions import Fraction as F
import mpmath as mp
mp.dps=530
r2=mp.sqrt(2)
# planner funcs
EPS=mp.mpf(10)**(-520)
def liP(k,x):
    if x==-1: return -(1-mp.mpf(2)**(1-k))*mp.zeta(k)
    tot,xn,n=mp.mpf(0),mp.mpf(1),1
    while True:
        xn*=x; tot+=xn/mp.mpf(n)**k
        if abs(xn)<EPS: return tot
        n+=1
def li32P(x,y):
    tot,H,xn,ym,n=mp.mpf(0),mp.mpf(0),mp.mpf(1),mp.mpf(1),1
    while True:
        xn*=x; tot+=xn*H/mp.mpf(n)**3
        if n>50 and abs(xn)<EPS: return tot
        ym*=y; H+=ym/mp.mpf(n)**2; n+=1
# my check_supplement funcs
def liM(k,x):
    N=int(450/float(-mp.log10(abs(x))))+60
    s,xp=mp.mpf(0),x
    for j in range(1,N+1):
        s+=xp/mp.mpf(j)**k; xp*=x
    return s
def li32M(x,y):
    N=int(450/float(-mp.log10(abs(x))))+60
    s,xn,yj,h=mp.mpf(0),x,mp.mpf(1),mp.mpf(0)
    for n in range(2,N+1):
        yj*=y; h+=yj/mp.mpf(n-1)**2; xn*=x
        s+=xn/mp.mpf(n)**3*h
    return s
# compare li32 on the four N5 args
for x,y in ((r2-1,1),(1/r2,1),(1-r2,1),(r2-1,-1)):
    d=li32P(x,y)-li32M(x,y)
    print('li32 x=%s y=%s P-M=%s'%(mp.nstr(x,4),mp.nstr(y,3),mp.nstr(d,3)))
# compare li on all 12 args
def fld(p,q): return mp.mpf(p.numerator)/p.denominator+r2*mp.mpf(q.numerator)/q.denominator
args=[(F(-1),F(1)),(F(1),F(-1)),(F(-3),F(2)),(F(0),F(1,2)),(F(1,2),F(0)),(F(1),F(-1,2)),
      (F(2),F(-1)),(F(-2),F(2)),(F(1,2),F(-1,4)),(F(1,2),F(1,4)),(F(1,2),F(-1,2)),(F(1,2),F(-3,8))]
for pq in args:
    x=fld(*pq)
    for k in (2,3,4,5):
        d=liP(k,x)-liM(k,x)
        if abs(d)>mp.mpf('1e-500'): print('li%d x=%s P-M=%s'%(k,mp.nstr(x,6),mp.nstr(d,3)))
print('li comparisons done')
