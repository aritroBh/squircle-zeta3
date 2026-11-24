from mpmath import mp, mpf, pi, sqrt
mp.dps = 50
wb = expv = None
j=1j
wb_ = mp.e**(j*pi/4); wb = wb_.conjugate()
worst = {}
for N in [400, 3000, 100000]:
    N = mpf(N); y = sqrt(pi/(2*N))
    for m in range(2, 8):
        lm = (-1)**(m+1) * y**m/m * ((N-1) - N/(1+y*wb)**m)
        Lm = N*y**(m+1) + y**m/m
        r = abs(lm)/Lm
        worst[m] = max(worst.get(m, 0), r)
for m in sorted(worst):
    print(f"m={m}: max |l_m|/L_m over N in {{400,3000,1e5}} = {mp.nstr(worst[m],10)}  (<=1 needed)")
