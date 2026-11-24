# step_00.py - D13h: verify algebraic identities for the final analytic cell
from flint import arb, acb, ctx
ctx.dps = 60
pi = arb.pi(); rt2 = arb(2).sqrt(); sqpi = pi.sqrt()
b = acb(0,1)*pi/4; b = b.exp()
I1 = acb(0,1)
a = pi/(2*rt2)

def Hc(s): return sqpi/4*(s/2).gamma()/(s/2+arb("1.5")).gamma()
def gF(x): return Hc(b*x)/x

def mufun(z):  # Binet mu via logGamma remainder
    return z.loggamma() - (z-arb("0.5"))*z.log() + z - (2*pi).log()/2

print("=== (A) Binet mu bound: |mu(z)| <= 1/(12|z|cos(arg z)) for arg=pi/4 ray ===")
worst = 0
for xv in [0.1, 0.3, 1.0, 3.0, 10.0, 100.0, 1e4, 1e6]:
    z = b*xv/2
    mu = mufun(z)
    th = arb(pi/4)
    bd = 1/(12*abs(z)*arb(th).cos())
    r = float(abs(mu)/bd)
    worst = max(worst, r)
    print(f"  x={xv:>9}: |mu|={float(abs(mu)):.3e}, bound={float(bd):.3e}, ratio={r:.4f}")
print("  worst ratio:", worst)

print("=== (B) g0 factored: |g(x)| <= (sqpi/4) 2^{3/2} x^{-5/2} e^{ReE}, ReE <= -t*rt2/24 + (459/24) t^3/(1-3t), t=1/x ===")
for xv in [50.0, 1e3, 1e5, 7.9788e4]:
    t = 1/xv
    s = b*xv
    u = s/2
    E = arb("1.5") - (u+1)*(1+3/(2*u)).log() + mufun(u) - mufun(u+arb("1.5"))
    ReE_up = -t*rt2/24 + (arb(459)/24)*t**3/(1-3*t)
    lhs = abs(gF(xv))
    rhs = sqpi/4*2**arb("1.5")*xv**arb("-2.5")*(ReE_up.exp())
    print(f"  x={xv:>10}: |g|={float(lhs):.6e}, ReE={float(acb(E).real):.3e} <= {float(ReE_up):.3e}? {acb(E).real <= ReE_up}, |g|/rhs={float(lhs/rhs):.5f}")

print("=== (C) z-factorization identity: z == e^{i*pi/4} rho zt3 e^E e^{mS} e^{om^2 invp2/2} e^{i phi} ===")
# replicate step_38 cell quantities at a point y (sigma=0.3), tiny box
y = arb("0.0627"); sig = arb("0.3")
rho = 1 + y*y*sig; t = rho*y; x0 = 1/t; N = pi/(2*y*y)
den1 = 1+rt2*rho*y+rho*rho*y*y
Omn = rt2*rho + rho*rho*y - 2*y*sig - y**3*sig**2
om = (pi/(2*rt2))*y*Omn/den1
numR = t**3 + 2*rt2*t**2 + 4*t + rt2
D2 = t**2 + 2*rt2*t + 4
invp2 = -(2/(pi*rho**3*y))*acb(numR, t+rt2)/D2
# step_38 z-factor z1*z3 with g0box -> exact gF(x0)
g0 = gF(x0)
z1z3 = (2*N*I1/pi**arb("1.5"))*g0*((-2*pi*invp2).sqrt())
# my factored form: e^{i pi/4} rho zt3 e^E  (zt3 = sqrt(Q3), Q3 = (numR+i(t+rt2))*2 e^{-i pi/4}/D2)
u = b*x0/2
E = arb("1.5") - (u+1)*(1+3/(2*u)).log() + mufun(u) - mufun(u+arb("1.5"))
Q3 = (acb(numR, t+rt2))*2*(-I1*pi/4).exp()/D2
zt3 = Q3.sqrt()
z1z3_mine = (I1*pi/4).exp()*rho*zt3*E.exp()
print("  z1z3      =", z1z3)
print("  z1z3 mine =", z1z3_mine)
print("  ratio:", z1z3/z1z3_mine)
