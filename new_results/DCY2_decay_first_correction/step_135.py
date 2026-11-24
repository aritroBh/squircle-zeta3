from mpmath import mp, mpf, pi, sqrt, exp, erf
mp.dps = 50
bR = pi/(2*sqrt(2))
for N in [400, 2000, 10000]:
    y = sqrt(pi/(2*N))
    c2 = 0.9*bR*y
    ReL1 = y*(1+bR*y)
    a = ReL1; c = c2
    tR = a/(2*c)
    efac = exp(a**2/(4*c))
    # exact full two-sided moment: 2 efac ∫_0^inf t^j e^{-c(t-tR)^2} dt
    for j in [0, 4, 6, 12]:
        exact = 2*efac*mp.quad(lambda t: t**j*exp(-c*(t-tR)**2), [0, mp.inf])
        # my Mtilde: efac * sum binom(j,i) tR^{j-i} c^{-(i+1)/2} Gamma((i+1)/2)
        from mpmath import gamma, binomial
        Mt = efac*sum(binomial(j,i)*tR**(j-i)*c**(-(i+1)/2)*gamma((i+1)/2) for i in range(j+1))
        # correction 2 efac tR^{j+1}/(j+1)
        corr = 2*efac*tR**(j+1)/(j+1)
        print(f"N={N} j={j}: exact={mp.nstr(exact,8)}, Mtilde={mp.nstr(Mt,8)}, M/exact={mp.nstr(Mt/exact,6)}, corr={mp.nstr(corr,6)}")
