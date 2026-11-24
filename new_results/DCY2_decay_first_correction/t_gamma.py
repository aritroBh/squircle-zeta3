
from flint import arb, ctx
ctx.dps = 40
PI = arb.pi()
def Gamma_upper_half(n, x):
    res = PI.sqrt()*(x.sqrt()).erfc()
    G = res
    for k in range(1, n+1):
        G = (k-arb(1)/2)*G + x**(k-arb(1)/2)*(-x).exp()
    return G
def Gamma_upper_int(n, x):
    tot = arb(0)
    for i in range(n):
        tot += x**i/arb(i).fac()
    return arb(n-1).fac()*(-x).exp()*tot
for s2, x in [(13, "0.742"), (13, "3.08"), (3, "1.0"), (9, "4.757"), (1, "0.5"), (13, "0.28"), (7, "12.0")]:
    n = (s2-1)//2
    a = arb(x)
    Gh = Gamma_upper_half(n, a)
    Gi = Gamma_upper_int(n, a) if n+1 == (s2+1)//2 and (s2 % 2 == 1) else None
    print(f"Gamma({s2}/2, {x}): half-int impl = {Gh.str(12)}")
for n, x in [(7, "0.742"), (4, "3.08"), (2, "1.0"), (1, "0.5")]:
    a = arb(x)
    print(f"Gamma({n}, {x}): int impl = {Gamma_upper_int(n, a).str(12)}")
