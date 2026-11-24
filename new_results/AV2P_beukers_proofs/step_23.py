from mpmath import mp, mpf, zeta, log, pi, nprint, binomial, polylog, clsin, pslq, catalan, sqrt
mp.dps = 210

K = 800
S1 = mpf(0); S1p = mpf(0); S2 = mpf(0); S3 = mpf(0)
Hk = mpf(0)
for k in range(0, K):
    if k > 0:
        Hk += 1/mpf(k)
    H2k = Hk + sum(1/mpf(k+j) for j in range(1, k+1)) if k>0 else mpf(0)
    c = binomial(2*k,k)/mpf(16)**k/(2*k+1)
    S1 += c*Hk; S1p += c*H2k; S2 += c/(2*k+1); S3 += c
print("S1 =", nprint(S1, 60))
print("S1' =", nprint(S1p, 60))
print("S2 - Cl2(pi/3) =", nprint(S2 - clsin(2,pi/3), 6))
print("S3 - pi/3 =", nprint(S3 - pi/3, 6))

odd = mpf(open('odd2.txt').read())
# odd part = (pi/8)(2 ln2 S3 + S1 - S1' + S2)
odd_check = pi/8*(2*log(2)*S3 + S1 - S1p + S2)
print("odd reconstruction residual:", nprint(odd - odd_check, 6))

# PSLQ basis
z3 = zeta(3); L2 = log(2); Cl = clsin(2, pi/3); G = catalan
Li3_half = polylog(3, mpf('0.5')); Li2_half = polylog(2, mpf('0.5'))
basis = [1, pi, pi**2, pi**3, z3, pi*z3, L2, pi*L2, pi**2*L2, L2**2, pi*L2**2, L2**3,
         Cl, pi*Cl, G, pi*G, Li3_half, Li2_half, pi*Li2_half, z3*L2, Cl*L2, log(3), pi*log(3),
         polylog(3, mpf('0.25')), pi**2*G if False else z3*pi**2]
names = ['1','pi','pi2','pi3','z3','pi*z3','log2','pi*log2','pi2*log2','log2^2','pi*log2^2',
         'log2^3','Cl2(pi/3)','pi*Cl2','G','pi*G','Li3(1/2)','Li2(1/2)','pi*Li2(1/2)',
         'z3*log2','Cl2*log2','log3','pi*log3','Li3(1/4)','pi2*z3']

def hunt(target, label, tol=None):
    v = [target] + basis
    rel = pslq(v, tol=tol or mpf('1e-190'), maxcoeff=mpf(10)**30, maxsteps=10000)
    print(f"--- PSLQ {label} ---")
    if rel is None:
        print("  no relation")
    else:
        print("  relation:", rel[0], "*target +", )
        for c,n in zip(rel[1:], names):
            if c: print(f"    {c} * {n}")
    return rel

r1 = hunt(odd, "odd part")
r2 = hunt(S1, "S1")
r3 = hunt(S1p, "S1'")
r4 = hunt(S1 - S1p, "S1-S1'")
with open('progress.log','a') as fh:
    fh.write("step_21: S1,S1',S2,S3 at 200d; S2=Cl2(pi/3),S3=pi/3 confirmed; PSLQ odd/S1/S1' no relation found?\n")
