from mpmath import mp, mpf, log, pi, nprint, quad, binomial, sqrt
mp.dps = 60
def G_series(y, K=500):
    y = mpf(y); S = mpf(0); Hk = mpf(0)
    for k in range(0, K):
        if k>0: Hk += 1/mpf(k)
        H2k = Hk + sum(1/mpf(k+j) for j in range(1, k+1)) if k>0 else mpf(0)
        S += binomial(2*k,k)*H2k*y**k
    return S
for y0 in ['0.02','0.05','0.0625']:
    y = mpf(y0); c = sqrt(1-4*y)
    Gs = G_series(y)
    c1 = (1/c)*log((1+c)/(2*c**2))
    c2 = (1/c)*log((1+c)/(2*c))
    print(f"y={y0}: series={nprint(Gs,25)}")
    print(f"   (1/c)log((1+c)/2c^2)={nprint(c1,25)}  diff={nprint(Gs-c1,4)}")
    print(f"   (1/c)log((1+c)/2c)  ={nprint(c2,25)}  diff={nprint(Gs-c2,4)}")
