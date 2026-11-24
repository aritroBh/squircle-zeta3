# step_13.py - validate final-cell W2/W3 segmented z-sums against pointwise TRUE integrals at N=1e9
from flint import arb, acb, ctx
ctx.dps = 60
pi = arb.pi(); rt2 = arb(2).sqrt(); sqpi = pi.sqrt()
def Hc(s): return sqpi/4*(s/2).gamma()/(s/2+arb("1.5")).gamma()
b = (acb(0,1)*pi/4).exp()
def gF(x): return Hc(b*x)/x
def Phi3b(x,N): return 2*N*(1+3*x+3*x**2)/(x**3*(1+rt2*x+x*x)**arb("1.5"))
def Phi4b(x,N): return 6*N*(1+4*x+6*x**2+4*x**3)/(x**4*(1+rt2*x+x*x)**2)
def G1bound(x):
    t = x/(2*rt2); return arb("0.75")*(1/t**2+1/t)+1/x
def G2bound(x):
    t = x/(2*rt2); return arb("0.375")*(2/t**3+1/t**2)+1/x**2
def h3(T): return rt2/3/(1+rt2*T+T*T)

N = arb(10)**9
y = (pi/(2*N)).sqrt(); sig = arb("0.3"); rho = 1+y*y*sig; t = rho*y; x0 = 1/t
numR = t**3+2*rt2*t**2+4*t+rt2
c2 = (pi/4)*rho**3*y*numR/(1+rt2*t+t*t)**2
c_w = arb(1); w = (N.log()+c_w).sqrt(); d = w/(2*c2).sqrt()
mS = (N*t**3/2)*h3(t)
P0t = abs(gF(x0))*mS.exp()
# TRUE W2 integral (both sides), sharp integrand, fine grid
M = 4000
W2t = arb(0)
for side in [1,-1]:
    du = d/M
    for j in range(M):
        u = (j+arb("0.5"))*du   # midpoint rule
        x = x0 + side*u
        integ = (-c2*u*u).exp()*(Phi4b(x,N)*u**4/24 + G2bound(x)*u**2/2)
        W2t += integ*du
W2t *= P0t
print("W2_true*N^1.5   =", float(W2t*N**arb("1.5")))
# my z-space sum value at this point (J=48, zstar=3, pointwise sharp)
c2reg = c2/y
zmax = d*rho*y/y**arb("0.5")  # = d rho sqrt(y)?? careful: u = z/(rho sqrt(y)) -> z = u rho sqrt(y)
print("zmax point =", float(d*rho*y**arb("0.5")))
# true integral in z form: W2 = y^3 rho^{3/2} P0reg ∫ e^{-K2reg z^2}(R4breg_z+Ga2breg_z) dz over z in [0, zmax]
P0reg = abs(gF(x0))*mS.exp()/t**arb("2.5")
K2reg = c2reg/rho**2
M2 = 4000
zmaxp = d*rho*y**arb("0.5")
S2true = arb(0)
for side in [1,-1]:
    dz = zmaxp/M2
    for j in range(M2):
        z = (j+arb("0.5"))*dz
        eta = 1 + side*z*y**arb("0.5")
        u = z/(rho*y**arb("0.5"))
        R4 = Phi4b(eta*x0, N)*u**4/24
        Ga2 = G2bound(eta*x0)*u**2/2
        S2true += (-K2reg*z*z).exp()*(R4+Ga2)*du*0 + 0  # placeholder
# simpler: directly verify W2 = y^3 rho^1.5 P0reg S2z with S2z from my sum -> compare W2/y^3
W2mine = (arb("9.3001")/(pi/2)**arb("1.5"))  # my W2reg approx from final cell (limit)
print("W2/y^3 point true =", float(W2t/y**3), " my limit W2reg = 4.724")
print("ratio true/mine =", float(W2t/y**3)/4.724)
