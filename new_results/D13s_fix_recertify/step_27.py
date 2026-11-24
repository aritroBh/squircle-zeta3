
import mpmath as mp
mp.mp.dps = 50
sqrt2 = mp.sqrt(2); pi = mp.pi

def Hfun(t):
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t))
    psi = sqrt2*t - mp.log(1+sqrt2*t+t*t)
    G = psi - (sqrt2/3)*t**3/(1+sqrt2*t+t*t)
    return G - sqrt2*t*(rho-1)**2/rho**2

# scan for sign of H over t in [1e-6, 1e6]
prev = None; prevt = None
crossings = []
ts = [mp.mpf(10)**(mp.mpf(k)/20) for k in range(-120, 241)]
for t in ts:
    h = Hfun(t)
    if prev is not None and prev*h < 0:
        crossings.append((prevt, t))
    prev, prevt = h, t
print("sign crossings of H between:", [(mp.nstr(a,4), mp.nstr(b,4)) for a,b in crossings])
for a,b in crossings:
    rt = mp.findroot(Hfun, (a+b)/2)
    print("zero of H at t* =", mp.nstr(rt, 12))
    # corresponding N on the saddle branch
    rho = mp.sqrt(1 + rt*rt/(1+sqrt2*rt)); y = rt/rho; N = pi/(2*y*y)
    print("   corresponding N =", mp.nstr(N, 12))
# values at key points
for tv in ['0.12633','0.127','0.13','0.5','1.0','10','1e3']:
    print(f"H({tv}) =", mp.nstr(Hfun(mp.mpf(tv)), 8))
# H/t^4 at small t
for tv in ['1e-2','1e-4','1e-8']:
    t = mp.mpf(tv)
    print(f"H({tv})/t^4 =", mp.nstr(Hfun(t)/t**4, 10))
