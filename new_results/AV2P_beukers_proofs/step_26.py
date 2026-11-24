from mpmath import mp, mpf, zeta, log, pi, nprint, catalan, clsin, binomial
mp.dps = 100
odd = mpf(open('odd2.txt').read())
K = 800
S1 = mpf(0); S1p = mpf(0); S2 = mpf(0); S3 = mpf(0); Hk = mpf(0)
for k in range(0, K):
    if k>0: Hk += 1/mpf(k)
    H2k = Hk + sum(1/mpf(k+j) for j in range(1, k+1)) if k>0 else mpf(0)
    c = binomial(2*k,k)/mpf(16)**k/(2*k+1)
    S1 += c*Hk; S1p += c*H2k; S2 += c/(2*k+1); S3 += c
Cl = clsin(2, pi/3)
print("8G/3 =", nprint(8*catalan/3, 40))
print("2 pi log2/3 + (S1-S1') + Cl =", nprint(2*pi*log(2)/3 + (S1-S1p) + Cl, 40))
print("odd =", nprint(odd, 40))
print("(pi/8)(2log2 S3 + S1 - S1' + S2) =", nprint(pi/8*(2*log(2)*S3 + S1 - S1p + S2), 40))
print("pi G/3 =", nprint(pi*catalan/3, 40))
