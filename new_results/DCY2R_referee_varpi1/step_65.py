from mpmath import mp, mpf, pi, sqrt
mp.dps = 50
j=1j
wb = mp.e**(-j*pi/4)
worst = {m: (0,None) for m in range(2,8)}
for i in range(2000):
    y = mpf("0.0627")*i/1999 + mpf("1e-9")
    for m in range(2,8):
        r = abs(1-(1+y*wb)**(-m))/(m*y)
        if r > worst[m][0]: worst[m] = (r, y)
for m in sorted(worst):
    print(f"m={m}: max |1-(1+y wb)^-m|/(m y) = {mp.nstr(worst[m][0],10)} at y={mp.nstr(worst[m][1],6)}")
# and |l_m|/L_m on a dense N grid for m=3,4,5
for m in [3,4,5]:
    w = 0
    for i in range(400, 20000, 37):
        N = mpf(i); y = sqrt(pi/(2*N))
        lm = (-1)**(m+1)*y**m/m*((N-1)-N/(1+y*wb)**m)
        Lm = N*y**(m+1)+y**m/m
        w = max(w, abs(lm)/Lm)
    print(f"m={m}: max |l_m|/L_m on dense grid = {mp.nstr(w,10)}")
