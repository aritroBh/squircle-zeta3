
import importlib, math, sys
import step_38 as S
from flint import arb
y = S.iv(math.sqrt(math.pi/(2*1.5e6)), math.sqrt(math.pi/(2*1e6)))
sig = S.iv(arb(0), arb("0.5"))
rho = 1 + y*y*sig
y0 = rho*y; x0 = 1/y0
x0lo = arb(x0.lower()); x0hi = arb(x0.upper())
N = S.pi/(2*y*y)
t = rho*y
numR = t**3 + 2*S.rt2*t**2 + 4*t + S.rt2
Rep2 = -(S.pi/2)*rho**3*y*numR/(1+S.rt2*t+t**2)**2
c2 = arb((-Rep2/2).lower()); sg = 1/(2*c2).sqrt()
w = (N.log()+arb("1.0")).sqrt(); d = w*sg
xr = x0hi + d
etaR = S.iv((xr/x0hi).lower(), (xr/x0lo).upper())
def show(name, e):
    print(f"{name}: mid={e.mid()} rad={e.rad()}")
show("etaR", etaR)
show("t", t)
tt = etaR*etaR + S.rt2*etaR*t + t*t
show("eta^2+..", tt)
den = S.rt2*etaR*tt
show("rt2*eta*(..)", den)
num = rho*(etaR + t*S.rt2)
show("rho*(eta+trt2)", num)
Gh = num/den
show("Gh", Gh)
lamR = S.a - (S.pi*rho/2)*Gh
show("lamR", lamR)
