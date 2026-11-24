
import importlib, sys, math
import step_38; importlib.reload(step_38)
from flint import arb
import step_38 as S

# instrument: replicate cell() up to lamR for [1e6, 1.5e6], sig sub-box 0
y = S.iv(S.ys if hasattr(S,'ys') else math.sqrt(math.pi/(2*1.5e6)), math.sqrt(math.pi/(2*1e6)))
sig = S.iv(arb(0), arb("0.5"))
rho = 1 + y*y*sig
y0 = rho*y; x0 = 1/y0
x0lo = arb(x0.lower()); x0hi = arb(x0.upper())
N = S.pi/(2*y*y)
c_w = arb("1.0")
t = rho*y
numR = t**3 + 2*S.rt2*t**2 + 4*t + S.rt2
Rep2 = -(S.pi/2)*rho**3*y*numR/(1+S.rt2*t+t**2)**2
c2 = arb((-Rep2/2).lower())
sg = 1/(2*c2).sqrt()
w = (N.log()+c_w).sqrt()
wsig = w*sg
d = wsig if arb(wsig.upper()) <= arb((x0lo/3).lower()) else x0lo/3
print("x0 box:", x0lo, x0hi, " rel width:", float((x0hi-x0lo)/x0lo))
print("d box:", d, " d/x0hi:", float((arb(d.upper())/x0lo)))
xr = x0hi + d
etaR = S.iv((xr/x0hi).lower(), (xr/x0lo).upper())
print("etaR:", etaR)
Ghv = S.Gh(etaR, rho, y)
lamR = S.a - (S.pi*rho/2)*Ghv
print("Gh box:", Ghv)
print("2a/(pi rho) =", float((2*S.a/(S.pi*rho)).mid()), " Gh(etaR.lo) should be < this")
print("lamR:", lamR)
print("y box:", y, "rho box:", rho)
