from mpmath import mp, mpf, pi, sqrt, log
mp.dps = 50
bR = pi/(2*sqrt(2)); SQ2 = sqrt(2)

def c2_true(N):
    N=mpf(N); y = sqrt(pi/(2*N))
    l2 = -(N-1)*y**2/2 + N*y**2/(2*(1+y*(1-1j)/SQ2)**2)
    return -l2.real, y

# 1. c2 >= 0.9 bR y on [400, 1e20]?
worst = 99
for N in [400, 1000, 1e4, 1e5, 1e6, 1e9, 1e12, 1e20]:
    c2, y = c2_true(N)
    r = c2/(bR*y)
    worst = min(worst, r)
    print(f"c2/(bR y) at N={N:10.0e}: {mp.nstr(r,12)}")
print("min:", mp.nstr(worst,8))

# 2. tail rates across the range
def check_rates(N):
    N=mpf(N); y = sqrt(pi/(2*N))
    kap2 = mpf('2.48')**2 + 4*log(N/400); kap = sqrt(kap2)
    c2m = mpf('0.9')*bR*y; sig = 1/sqrt(2*c2m); D = kap*sig
    w = y*D; u = 1+w
    Q = u**2+SQ2*y*u+y**2
    # true -Lambda'(Delta)
    LamP = y*((N-1)/u - N*(2*u+SQ2*y)/(2*Q)) - bR
    # certified rDelta
    rD = bR*y*(2*D+y*D**2-mpf('0.514'))/(1+w)**2
    # eLamD bound: Lambda(Delta) <= -kap^2/4 + kap*1.117*sqrt(y)
    Lam = (N-1)*log(u) + (N/2)*log(1+SQ2*y+y**2) - (N/2)*log(Q) - bR*D
    ebound = -kap2/4 + kap*mpf('1.117')*sqrt(y)
    # left rate: min over u in (0,uD] of P_true(u)/(2uQ) vs lamL_code
    uD = 1-w
    def P_true(uu):
        Qq = uu**2+SQ2*y*uu+y**2
        return 2*(N-1)*y*Qq - N*y*uu*(2*uu+SQ2*y) - 2*bR*uu*Qq
    def P_code(uu):
        return -2*bR*uu**3 - (pi+2)*y*uu**2 + (SQ2*(N-2)*y**2 - 2*bR*y**3)*uu + 2*(N-1)*y**3
    # check P_code == P_true?
    diff = max(abs(P_true(mpf(k)/20)-P_code(mpf(k)/20)) for k in range(21))
    QD = uD**2+SQ2*y*uD+y**2
    lamL_code = min(P_code(uD), pi*y/2)/(2*uD*QD)
    # true min rate on (0,uD]
    import mpmath
    grid = [uD*k/200 for k in range(1,201)]
    truerate = min(P_true(uu)/(2*uu*(uu**2+SQ2*y*uu+y**2)) for uu in grid)
    return dict(LamP=LamP, rD=rD, Lam=Lam, ebound=ebound, lamL=lamL_code, truerate=truerate, Pdiff=diff, w=w)

print()
for N in [400, 3000, 1e4, 1e5, 1e6]:
    r = check_rates(N)
    print(f"N={N:8.0f}: -Lam'(D)={mp.nstr(r['LamP'],6)} rDelta={mp.nstr(r['rD'],6)} ok={r['rD']<=r['LamP'] and r['rD']<=bR} | Lam(D)={mp.nstr(r['Lam'],7)} <= {mp.nstr(r['ebound'],7)} ok={r['Lam']<=r['ebound']}")
    print(f"          lamL_code={mp.nstr(r['lamL'],6)} truerate={mp.nstr(r['truerate'],6)} ok={r['lamL']<=r['truerate']}  |P_true-P_code|={mp.nstr(r['Pdiff'],4)}")
