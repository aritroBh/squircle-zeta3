from mpmath import mp, mpf, zeta, log, pi, nprint, quad, clsin, binomial, sqrt, pslq, catalan
mp.dps = 60
def g(y, K=800):
    y = mpf(y); S = mpf(0)
    for n in range(1, K):
        H2n = sum(1/mpf(k) for k in range(1,2*n+1))
        Hn = sum(1/mpf(k) for k in range(1,n+1))
        S += (H2n-Hn)*y**n/(n*binomial(2*n,n))
    return S
for y0 in ['0.25','0.09','0.16']:
    gy = g(y0)
    print(f"g({y0}) =", nprint(gy, 30))
    y = mpf(y0)
    # candidate forms: 2 asin^2(sqrt y/2)?  f? something with log and asin
    cands = {
        '4 asin^2(sqrt(y)/2)': 4*mp.asin(sqrt(y)/2)**2,
        '2 asin^2(sqrt y /2)': 2*mp.asin(sqrt(y)/2)**2,
        'asin^2(sqrt y)': mp.asin(sqrt(y))**2,
    }
    for nm, v in cands.items():
        if abs(gy-v) < mpf('1e-40'):
            print("   MATCH:", nm)
# numeric integral check: g(y) = int_0^1 [f0(xy)-f0(x^2 y)]/(1-x) dx
def f0(z):
    z = mpf(z)
    return 2*sqrt(z)*mp.asin(sqrt(z)/2)/sqrt(4-z)
for y0 in ['0.25','0.09']:
    y = mpf(y0)
    gi = quad(lambda x: (f0(x*y)-f0(x*x*y))/(1-x), [0,1])
    print(f"g({y0}) vs integral:", nprint(g(y0)-gi, 4))
# PSLQ for g(1/4) against arcsin/log stuff
G = catalan
for y0 in ['0.25']:
    y = mpf(y0); gy = g(y)
    a = mp.asin(sqrt(y)/2)
    basis = [1, a, a**2, log(2-sqrt(4-y)) if (4-y)>0 else 0, log(2), log(3), pi, pi**2]
    r = pslq([gy]+basis, tol=mpf('1e-50'), maxcoeff=10**10)
    print("g(1/4) pslq:", r)
# Try: g(y) = 4 asin(sqrt y/2)^2 - something?
y = mpf('0.25'); gy = g(y)
print("g - 4asin^2:", nprint(gy - 4*mp.asin(sqrt(y)/2)**2, 5))
print("g - 2asin^2:", nprint(gy - 2*mp.asin(sqrt(y)/2)**2, 5))
print("g/asin^2:", nprint(gy/(2*mp.asin(sqrt(y)/2)**2), 10))
