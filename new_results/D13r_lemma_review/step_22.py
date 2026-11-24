import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
def RePhi1(x, N): return mp.re((-1+1j)*a + N*(1/x - bb/(1+bb*x)))
for yv in ['0.32519','0.1','0.0789','0.03']:
    y = mp.mpf(yv); N = pi/(2*y*y); th = pi/(sqrt2*y)
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2)
    x0 = 1/(rho*y); t = rho*y
    S = mp.log(1+sqrt2*t+t*t) - sqrt2*t
    lhs = mp.re(Phi(x0,N)) + th
    rhs = -(N/2)*S
    print(f"y={yv}: N={mp.nstr(N,6)} rho={mp.nstr(rho,12)} RePhi'(x0)={mp.nstr(RePhi1(x0,N),6)}")
    print(f"   lhs={mp.nstr(lhs,16)}  rhs={mp.nstr(rhs,16)}  reldiff={mp.nstr(abs(lhs-rhs)/abs(lhs),6)}")
# Direct: numerically maximize RePhi to find the TRUE x0 independently (no saddle equation)
y = mp.mpf('0.32519'); N = pi/(2*y*y); th = pi/(sqrt2*y)
g = lambda x: mp.re(Phi(x,N))
x0num = mp.findroot(lambda x: RePhi1(x,N), 1/(y*(1+y*y/2)))
print("independent x0 from RePhi'=0:", mp.nstr(x0num,14), "vs 1/(rho y) =", mp.nstr(1/(mp.findroot(lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2)-1, 1+y*y/2)*y),14))
print("lhs at independent x0:", mp.nstr(mp.re(Phi(x0num,N))+th,16))
