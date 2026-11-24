from flint import arb, acb, ctx
ctx.dps = 60
pi = arb.pi(); rt2 = arb(2).sqrt(); sqpi = pi.sqrt()
b = (acb(0,1)*pi/4).exp()
I1 = acb(0,1)
a = pi/(2*rt2)

def Hc(s): return sqpi/4*(s/2).gamma()/(s/2+arb("1.5")).gamma()
def gF(x): return Hc(b*x)/x
def mufun(z): return z.lgamma() - (z-arb("0.5"))*z.log() + z - (2*pi).log()/2

print("=== (A) Binet mu bound: |mu(z)| <= 1/(12|z|cos(pi/4)) on ray arg z = pi/4 ===")
worst = 0
for xv in [0.1, 0.3, 1.0, 3.0, 10.0, 100.0, 1e4, 1e6]:
    z = b*xv/2
    mu = mufun(z)
    bd = 1/(12*abs(z)*arb(pi/4).cos())
    r = float(abs(mu)/bd); worst = max(worst, r)
    print(f"  x={xv:>9}: |mu|={float(abs(mu)):.3e}, bound={float(bd):.3e}, ratio={r:.4f}")
print("  worst ratio:", worst)

print("=== (B) g0 factored bound ===")
for xv in [50.0, 1e3, 1e5, 7.9788e4]:
    t = 1/xv
    s = b*xv; u = s/2
    E = arb("1.5") - (u+1)*(1+3/(2*u)).log() + mufun(u) - mufun(u+arb("1.5"))
    ReE = acb(E).real
    ReE_up = -t*rt2/24 + (arb(459)/24)*t**3/(1-3*t)
    lhs = abs(gF(xv))
    rhs = sqpi/4*2**arb("1.5")*xv**arb("-2.5")*ReE_up.exp()
    print(f"  x={xv:>10}: ReE={float(ReE):.6e} <= {float(ReE_up):.6e}? {ReE <= ReE_up}; |g|/rhs = {float(lhs/rhs):.6f}")

print("=== (C) z-factorization identity ===")
y = arb("0.0627"); sig = arb("0.3")
rho = 1 + y*y*sig; t = rho*y; x0 = 1/t; N = pi/(2*y*y)
numR = t**3 + 2*rt2*t**2 + 4*t + rt2
D2 = t**2 + 2*rt2*t + 4
invp2 = -(2/(pi*rho**3*y))*acb(numR, t+rt2)/D2
g0 = gF(x0)
z1z3 = (2*N*I1/pi**arb("1.5"))*g0*((-2*pi*invp2).sqrt())
u = b*x0/2
E = arb("1.5") - (u+1)*(1+3/(2*u)).log() + mufun(u) - mufun(u+arb("1.5"))
Q3 = (acb(numR, t+rt2))*2*(-I1*pi/4).exp()/D2
z1z3_mine = (I1*pi/4).exp()*rho*Q3.sqrt()*E.exp()
print("  z1z3      =", z1z3)
print("  z1z3 mine =", z1z3_mine)
print("  |ratio|:", abs(z1z3/z1z3_mine))

print("=== (D) check logcf / h(w) expansion constants: |3/2-(u+1)log(1+3/(2u))| <= 3t/4 + (459/24)t^3/(1-3t) ===")
for xv in [10.0, 100.0, 1e4, 8e4]:
    t = 1/xv; u = b*xv/2
    val = abs(arb("1.5") - (u+1)*(1+3/(2*u)).log())
    bd = 3*t/4 + (arb(459)/24)*t**3/(1-3*t)
    print(f"  x={xv:>8}: |val|={float(val):.6e} <= {float(bd):.6e}? {val <= bd}, ratio={float(val/bd):.4f}")
