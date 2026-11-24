import mpmath as mp
mp.mp.dps = 50

def f_a(a, t):
    if a == 0: return 1/t
    if a == 1: return 1/(1-t)
    return 1/(1+t)

def H1(a, t):
    if a == 0: return mp.log(t)
    if a == 1: return -mp.log(1-t)
    return mp.log(1+t)

def H2(a, b, x):
    return mp.quad(lambda t: f_a(a, t)*H1(b, t), [0, x])

lam = mp.log(2)
def C_1m1(x):   # candidate for H_{1,-1}
    return -mp.polylog(2, (1+x)/2) - mp.log(1-x)*mp.log(1+x) + lam*mp.log(1+x) - lam**2/2 + mp.pi**2/12
def C_m11(x):   # candidate for H_{-1,1}
    return mp.polylog(2, (1+x)/2) - lam*mp.log(1+x) + lam**2/2 - mp.pi**2/12

for xv in [mp.mpf('0.2'), mp.mpf('0.5'), mp.mpf('0.8'), mp.mpf('0.95')]:
    d1 = H2(1,-1,xv) - C_1m1(xv)
    d2 = H2(-1,1,xv) - C_m11(xv)
    print(f"x={xv}:  H(1,-1)-C={mp.nstr(d1,5)}   H(-1,1)-C={mp.nstr(d2,5)}")
