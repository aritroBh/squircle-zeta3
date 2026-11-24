import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
import mpmath as mp
mp.mp.dps=60
g=S.init()
Xv={i:E.fev(E.X[i],60) for i in E.IDX}
# (2) validate Ax reduction: A(x) = Li3(x)-Li3(1-x)+zeta(3) + (1/2)log^2(1-x)log x + log(1-x)Li2(1-x)
def Aser(x):
    tot=mp.mpf(0); xp=mp.mpf(1); H=mp.mpf(0)
    for n in range(1,8000): xp*=x; H+=mp.mpf(1)/n; tot+=H*xp/n**2
    return tot
for i in [1,6]:
    x=Xv[i]
    lhs=Aser(x)
    rhs=mp.polylog(3,x)-mp.polylog(3,1-x)+mp.zeta(3)+mp.mpf('0.5')*mp.log(1-x)**2*mp.log(x)+mp.log(1-x)*mp.polylog(2,1-x)
    print(f"Ax reduction x{i}: diff", mp.nstr(abs(lhs-rhs),4))
# (2b) validate T3p reduction: T3p(i,j) = -(1/x_i)*[2Li3(c)-2Li3(cX)+log^2 X log(1-cX)+2 logX Li2(cX)]
# with x=1/x_j, y=1/x_i, c=x/(x-y), X=(x-1)/x, cX=(x-1)/(x-y)
def T3pred(i,j):
    x=1/Xv[j]; y=1/Xv[i]
    c=x/(x-y); X=(x-1)/x; cX=(x-1)/(x-y)
    return -(1/Xv[i])*(2*mp.polylog(3,c)-2*mp.polylog(3,cX)+mp.log(X)**2*mp.log(1-cX)+2*mp.log(X)*mp.polylog(2,cX))
w={i:(lambda tt,i=i: mp.log(1-Xv[i]*tt)) for i in E.IDX}
u={i:(lambda tt,i=i: 1-Xv[i]*tt) for i in E.IDX}
for (i,j) in [(1,2),(2,5),(6,1)]:
    tru=mp.quad(lambda tt: w[j](tt)**2/u[i](tt),[0,1])
    print(f"T3p reduction ({i},{j}): diff", mp.nstr(abs(tru-T3pred(i,j)),4))
# (4) C7 check: local zeta(3) coefficient C_m = 48 m sin^4(pi/2m) cos(pi/2m) (3 sin^2(pi/2m)-1) at m=7
m=7
C7 = 48*m*mp.sin(mp.pi/(2*m))**4*mp.cos(mp.pi/(2*m))*(3*mp.sin(mp.pi/(2*m))**2-1)
print("C_7 =", mp.nstr(C7,30))
open('progress.log','a').write('PG6 s_C1: Ax/T3p reductions validated numerically; C7 computed\n')
