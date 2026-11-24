from mpmath import mp, mpf, pi, sqrt, log, exp
mp.dps = 50
bR = pi/(2*sqrt(2)); SQ2 = sqrt(2)

def check_tails(N):
    N=mpf(N); y = sqrt(pi/(2*N)); rho=1/y
    kap2 = mpf('2.48')**2 + 4*log(N/400); kap = sqrt(kap2)
    c2m = mpf('0.9')*bR*y; sig = 1/sqrt(2*c2m); D = kap*sig
    w = y*D
    def Lam(t):
        u = 1+y*t
        Q = u**2+SQ2*y*u+y**2
        return (N-1)*log(u) + (N/2)*log(1+SQ2*y+y**2) - (N/2)*log(Q) - bR*t
    def LamP(t):
        u = 1+y*t
        Q = u**2+SQ2*y*u+y**2
        return y*((N-1)/u - N*(2*u+SQ2*y)/(2*Q)) - bR
    rD = bR*y*(2*D+y*D**2-mpf('0.514'))/(1+w)**2
    # (a) rD <= min(-Lam'(D), bR)
    a_ok = (rD <= -LamP(D)) and (rD <= bR)
    # (b) rD <= -Lam'(t) for all t in [D, 30*rho] (grid) ; asymptotic for t beyond
    ts = [D + (30*rho-D)*k/400 for k in range(401)]
    b_ok = all(rD <= -LamP(t) for t in ts)
    # asymptotic: -Lam'(t) >= bR - pi/(sqrt2 u) - pi y/(2 u^2), u >= 30 -> >= rD ?
    u30 = 1+y*30*rho
    asym = bR - pi/(SQ2*30) - pi*y/(2*30**2)
    c_ok = (asym >= rD)
    # (c) left tail: Lambda(t) <= Lambda(-D) - lamL*(-D-t) for t in (-rho+1e-6, -D]
    uD = 1-w
    QD = uD**2+SQ2*y*uD+y**2
    P_c = -2*bR*uD**3 - (pi+2)*y*uD**2 + (SQ2*(N-2)*y**2 - 2*bR*y**3)*uD + 2*(N-1)*y**3
    lamL = min(P_c, pi*y/2)/(2*uD*QD)
    ts2 = [-D - (rho-1-D)*k/300 for k in range(301)]
    d_ok = all(Lam(t) <= Lam(-D) - lamL*(-D-t) + 1e-12 for t in ts2)
    # (d) eLamD: Lambda(D) <= -kap^2/4 + kap*1.117 sqrt(y)
    e_ok = Lam(D) <= -kap2/4 + kap*mpf('1.117')*sqrt(y)
    # (e) window-tail condition a <= c2 Delta/2 (for TailMom validity): |Re l1| <= c2 Delta/2
    ReL1 = y*(1+bR*y/(1+SQ2*y+y**2))
    l2re = -((N-1)*y**2/2 - N*y**2*(1+SQ2*y)/(2*(1+SQ2*y+y**2)**2))
    f_ok = ReL1 <= -l2re*D/2  # note -l2re = c2 = -Re l2 ... check sign
    return a_ok, b_ok, c_ok, d_ok, e_ok, f_ok, rD, lamL

for N in [400, 1000, 3000, 10000, 100000, 1000000]:
    a,b,c,d,e,f,rD,lamL = check_tails(N)
    print(f"N={N:8d}: rD<=min(-Lam'(D),bR):{a}  rD valid on [D,30rho]:{b}  asym:{c}  lefttail:{d}  eLamD:{e}  ReL1<=c2D/2:{f}")
