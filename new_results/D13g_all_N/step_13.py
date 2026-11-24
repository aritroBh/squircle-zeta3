
from flint import arb, acb, ctx
ctx.dps = 60
pi = arb.pi(); rt2 = arb(2).sqrt(); a = pi/(2*rt2)
b = (acb(0,1)*pi/4).exp(); I1 = acb(0,1)

def ReX2_old(x2, N, th):
    vf = rt2/x2 + 1/x2**2
    return -a*x2 + th - (N/2)*(1+vf).log()

def ReX2_reg(y, rho, eta):
    # = -a x + th - (N/2) log(1+rt2/x+1/x^2), x = eta/(rho y), N = pi/(2 y^2), th = pi/(rt2 y)
    u = rt2*rho*y/eta + (rho*y/eta)**2
    return (pi/(rt2*y))*(1 - eta/(2*rho)) - (pi/(4*y**2))*(1+u).log()

for Nv in [100, 400, 10**6]:
    N = arb(Nv)
    y = (pi/(2*N)).sqrt(); th = (pi*N).sqrt()
    # crude x0: solve ReP1=0 by bisection
    def ReP1(x): return -a + N/x - N*(x+1/rt2)/(1+rt2*x+x*x)
    lo, hi = arb(1), arb(10)**6
    for _ in range(300):
        mid=(lo+hi)/2
        if ReP1(mid)>0: lo=mid
        else: hi=mid
    x0=(lo+hi)/2; x2=x0/2
    rho = 1/(x0*y); eta = arb("0.5")
    v1 = ReX2_old(x2, N, th); v2 = ReX2_reg(y, rho, eta)
    print(f"N={Nv}: old={v1.mid()}  reg={v2.mid()}  diff={abs(v1-v2)}")
