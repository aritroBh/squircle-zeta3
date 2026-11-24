import mpmath as mp, time, json
mp.mp.dps = 50
t0 = time.time()

def S12(z):
    return mp.zeta(3) - mp.polylog(3, 1-z) + mp.log(1-z)*mp.polylog(2, 1-z) \
           + mp.mpf('0.5')*mp.log(1-z)**2*mp.log(z)
def cumulants(th):
    z = mp.sin(th/2)**2; w = 1/mp.tan(th/2)**2; ell = mp.log(mp.cos(th/2)**2); u = w*ell
    return (1+u, 1+mp.pi**2/3-2*mp.polylog(2,z)-u**2, -12*w*S12(z)+2*(1+u**3))

beta = mp.pi**2/8; m = mp.mpf(2)
data = json.load(open('contact_results.json'))
out = {}
for key, v in data.items():
    out[key] = {kk: mp.mpf(vv) for kk, vv in v.items()}
# keys are "('pi/3', 100000)" style
def get(thn, n): return out[str((thn, n))]

print("== leading-law ratio vs predicted correction 1 + ((m-1)log Lam + m k1)/(m^2 Lam) ==")
for thn, th in [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3)]:
    k1 = cumulants(th)[0]
    for n in [10**4, 10**5, 10**6, 10**7, 10**8]:
        d = get(thn, n)
        Nn = mp.sqrt(beta*mp.mpf(n)); Lam = mp.log(m*Nn)
        ratio = d['D']/d['L0']
        pred = 1 + ((m-1)*mp.log(Lam) + m*k1)/(m**2*Lam)
        rem = (ratio - pred)*(Lam/mp.log(Lam))**2   # O((log Lam/Lam)^2) remainder, rescaled
        print(f" th={thn:6s} n={n:8d}: D/L0={mp.nstr(ratio,7)} pred={mp.nstr(pred,7)} "
              f"(D/L0-pred)*(Lam/logLam)^2={mp.nstr(rem,6)}")

print("\n== full expansion in terms of log n (paper-friendly) ==")
# ratio = 1 + [(m-1) log log n + m k1 + log(m^m beta) - (m-1) log m]/(m log n) + O((loglog/log)^2)
for thn, th in [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3)]:
    k1 = cumulants(th)[0]
    for n in [10**5, 10**7, 10**8]:
        d = get(thn, n)
        L = mp.log(n)
        pred = 1 + ((m-1)*mp.log(L) + m*k1 + mp.log(m**m*beta) - (m-1)*mp.log(m))/(m*L)
        print(f" th={thn:6s} n={n:8d}: D/L0={mp.nstr(d['D']/d['L0'],7)} pred(logn form)={mp.nstr(pred,7)}")

print("\n== K(c) identity: Klayer(c) - (1-c^2)/4  vs  (1-c) k1(theta)/m ==")
def hc(s, c):
    if s < 1: return 1 + s - mp.sqrt(1 + 2*c*s + s**2)
    r = mp.sqrt(1 + 2*c*s + s**2)
    return (1-c) - (1-c**2)/(r + s + c)
def rho_st(s, c):
    r = mp.sqrt(1 + 2*c*s + s**2)
    return (1-c**2)/(r + s + c)
for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3]:
    c = mp.cos(th); k1 = cumulants(th)[0]
    i1 = mp.quad(lambda s: hc(s,c)/s, [0, 1])
    i2 = mp.quad(lambda s: rho_st(s,c)/s, [1, mp.inf])
    Klayer = (i1-i2)/2
    lhs = Klayer - (1-c**2)/4
    rhs = (1-c)*k1/m
    print(f" th={mp.nstr(th,6)}: Klayer-(1-c^2)/4={mp.nstr(lhs,12)}  (1-c)k1/m={mp.nstr(rhs,12)}  diff={mp.nstr(lhs-rhs,2)}")
print("elapsed", time.time()-t0)
