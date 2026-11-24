"""DCY2 error-budget prototype (mpmath, not certified). Every bound here is meant
to be replaced by an arb-certified one in verify_dcy2.py."""
from mpmath import mp, mpf, pi, sqrt, exp, gamma, log, factorial

j = 1j
omega = exp(j*pi/4); wb = exp(-j*pi/4)
bR = pi/(2*sqrt(2))          # Re b, b = pi*wb/2
b  = pi*wb/2

class P:  # parameters
    N1   = 400        # analytic range starts here
    kap0 = 2.5        # kappa^2 = kap0^2 + 4 log(N/N1)
    C2   = mpf(10)**6

def quant(N):
    N = mpf(N)
    q = {}
    q['N'] = N
    q['y'] = y = sqrt(pi/(2*N))
    q['rho'] = rho = 1/y
    q['theta'] = th = sqrt(pi*N)
    q['kap'] = kap = sqrt(P.kap0**2 + 4*log(N/P.N1))
    # ell_m exact
    def ell(m):
        if m == 1: return (N-1)*y - N*y/(1+y*wb) + j*pi*omega/2
        return (-1)**(m+1) * y**m/m * ((N-1) - N/(1+y*wb)**m)
    l = {m: ell(m) for m in range(1,7)}
    q['l1'], q['l2'] = l[1], l[2]
    # crude |ell_m| bounds for m >= 2
    Lb = lambda m: N*y**(m+1) + y**m/m
    # c2 = -Re l2 > 0 lower bound (certified separately); here: 0.9*bR*y
    q['c2'] = c2 = -l[2].real
    q['c2min'] = c2m = 0.9*bR*y
    assert c2 >= c2m, (N, c2, c2m)
    # sigma from c2min
    q['sig'] = sig = 1/sqrt(2*c2m)
    q['Delta'] = Delta = kap*sig
    q['yDelta'] = y*Delta
    assert q['yDelta'] < 0.75
    # S(y) sup
    S_sup = (pi*y/6)/(1-y)
    q['eRe_phi'] = y*exp(-th)*exp(S_sup)          # upper bound of |e^{phi(s0)}| real factor... e^{Re phi} <= y e^{-th} e^{S_sup}
    # |H(s0)| bounds: |g| <= 1 + 3y/4 + Cg y^2 ; Cg certified separately (placeholder 1.0)
    Cg = mpf('1.0')
    g_up = 1 + 3*y/4 + Cg*y**2
    g_lo = 1 - 3*y/4 - Cg*y**2
    q['H_up'] = sqrt(pi/2)*rho**(-mpf(3)/2)*g_up
    q['H_lo'] = sqrt(pi/2)*rho**(-mpf(3)/2)*g_lo
    # |h_j| bounds (Watson; placeholder factors 1+y)
    H1 = mpf(3)/2*y*(1+y); H2 = mpf(15)/8*y**2*(1+y); H3 = mpf(35)/16*y**3*(1+y)
    # |ell_m| for m=3,4,5 exact values available; use exact (they're computed)
    L3, L4, L5 = abs(l[3]), abs(l[4]), abs(l[5])
    # ----- window series sups
    yD = q['yDelta']
    B3 = (N*y**4 + y**3/3)/(1-yD)              # |R3(t)| <= B3 |t|^3
    B3sup = B3*Delta**3
    B6 = (N*y**7 + y**6/6)/(1-yD)              # |R6(t)| <= B6 |t|^6
    B6sup = B6*Delta**6
    alpha = L3 + L4*Delta + L5*Delta**2        # |A(t)| <= alpha |t|^3
    beta2 = L3*L5 + L4**2/2 + L4*L5*Delta + L5**2*Delta**2/2   # |t|^8 coeff
    beta3 = (3*L3**2*L4 + (3*L3**2*L5+3*L3*L4**2)*Delta + (6*L3*L4*L5+L4**3)*Delta**2
             + (3*L3*L5**2+3*L4**2*L5)*Delta**3 + 3*L4*L5**2*Delta**4 + L5**3*Delta**5)/6
    # ----- M4 (H fourth-moment remainder) : |F - sum_0^3 h_j t^j| <= M4 |t|^4 /24
    # numerator int u^4 sqrt(2u) e^{-x u} du = sqrt2 Gamma(11/2) x^{-11/2}, x = 1/(sqrt2 y)
    x = (1/y-Delta)/sqrt(2)   # t<0 rate: e^{-u(1/y-Delta)/sqrt2}
    M4num = sqrt(2)*gamma(mpf(11)/2)*x**(-mpf(11)/2)
    M4 = M4num / q['H_lo']      # ~ C y^4 (1-yDelta)^{-11/2}
    # ----- |cdot|-moment bounds  Mtilde_j
    ReL1 = abs(l[1].real)
    tR = ReL1/(2*c2m)
    efac = exp(ReL1**2/(4*c2m))
    def Mtilde(jj):
        tot = 0
        for i in range(jj+1):
            tot += mp.binomial(jj,i)*tR**(jj-i)*c2m**(-(i+1)/2)*gamma((i+1)/2)
        return efac*tot
    # ----- E_mom coefficients: E(t) bounded by sum Bj |t|^j (times e^{B3|t|^3} factor)
    eB3 = exp(B3sup)
    eA  = exp(alpha*Delta**3)
    eB6 = exp(B6sup)
    Hsup = 1 + H1*Delta + H2*Delta**2 + H3*Delta**3
    # E_a : (M4 |t|^4/24) e^{B3|t|^3}
    # E_b : |1+h1t+h2t^2+h3t^3| * [beta2|t|^8 + beta3|t|^10 + alpha^4 e^{alpha Delta^3}|t|^12/24
    #        + e^{alpha Delta^3} B6 e^{B6sup}|t|^6] e^{B3|t|^3}   (Q used e^{R3} factor too)
    # E_drop : P_drop coefficients
    E_drop = {
        6:  H1*L5 + H2*L4 + H3*L3,
        7:  H2*L5 + H3*L4,
        8:  H1*L3*L4 + H2*L3**2/2 + H3*L5,
        9:  H2*L3*L4 + H3*L3**2/2,
        10: H1*L3**3/6,
        11: H2*L3**3/6,
        12: H3*L3**3/6,
    }
    V3 = max(mpf(0), B3*Delta**3 - c2m*Delta**2)
    VA = max(mpf(0), alpha*Delta**3 - c2m*Delta**2)
    eV3 = exp(V3); eVA = exp(VA)
    E_mom = (M4/24)*Mtilde(4)*eV3
    E_mom += Hsup*( beta2*Mtilde(8) + beta3*Mtilde(10) + alpha**4*eVA/24*Mtilde(12)
                    + eVA*B6*eB6*Mtilde(6) )
    E_mom += sum(cf*Mtilde(p) for p,cf in E_drop.items())
    q['E_mom'] = E_mom
    q['V3'] = V3; q['VA'] = VA
    # ----- window-tail moments TP_p
    # |p_p| bounds
    pp = {
        0: mpf(1),
        1: H1,
        2: H2,
        3: H3 + L3,
        4: abs(l[4]) + H1*L3,
        5: L5 + H1*L4 + H2*L3,
        6: L3**2/2,
        7: abs(l[3]*l[4]) + H1*L3**2/2,
        9: L3**3/6,
    }
    from mpmath import gammainc
    assert ReL1 <= c2m*Delta/2, (N, ReL1, c2m*Delta/2)
    def TailMom(p):
        # int_Delta^inf t^p e^{-c2 t^2 + |Re l1| t} dt <= int_Delta^inf t^p e^{-c2 t^2/2} dt
        # = (1/2)(c2m/2)^{-(p+1)/2} Gamma((p+1)/2, c2m Delta^2/2)
        return mpf(1)/2*(c2m/2)**(-(p+1)/2)*gammainc((p+1)/2, c2m*Delta**2/2, regularized=False)
    E_tailmom = sum(cf*TailMom(p) for p,cf in pp.items())
    q['E_tailmom'] = E_tailmom
    # ----- outer tails
    # M(Delta) and M(-Delta): log M = Lambda(t)
    def LogM(t):
        u = 1+y*t
        Q = u**2 + sqrt(2)*y*u + y**2
        return (N-1)*log(u) + (N/2)*log(1+sqrt(2)*y+y**2) - (N/2)*log(Q) - bR*t
    LM_D = LogM(Delta)
    LM_mD = LogM(-Delta)
    # right tail rate: M(t) <= M(Delta) e^{-rDelta (t-Delta)}, rDelta = min(-Lambda'(Delta), b)
    yD2 = y*Delta
    uD_R = 1+yD2
    QR = uD_R**2 + sqrt(2)*y*uD_R + y**2
    LamPR = y*((N-1)/uD_R - N*(uD_R+y/sqrt(2))/QR) - bR
    assert LamPR < 0
    rDelta = min(-LamPR, bR)
    lamR = rDelta
    # left rate lambda_L: min of P(u)/(2 u Q(u)) on u in (0, 1-yD] ; use value at u = 1-yD (certified separately)
    def Pofu(u):
        return -2*bR*u**3 - (pi+2)*y*u**2 + (sqrt(2)*(N-2)*y**2 - 2*bR*y**3)*u + 2*(N-1)*y**3
    uD = 1-yD
    QD = uD**2 + sqrt(2)*y*uD + y**2
    lamL = Pofu(uD)/(2*uD*QD)
    assert lamL > 0
    # |H| bounds on tails: right: sqrt(pi/2)((rho+Delta)/sqrt2)^{-3/2}; left: sqrt(pi/2)*2^{3/4} for x>=1
    CH_R = sqrt(pi/2)*((rho+Delta)/sqrt(2))**(-mpf(3)/2)
    CH_L = sqrt(pi/2)*2**(mpf(3)/4)
    T_R = q['eRe_phi']*exp(LM_D)*CH_R/lamR
    # left: [rho-1-Delta piece] + [0,1] piece
    T_L1 = q['eRe_phi']*exp(LM_mD)*CH_L/lamL
    # [0,1] piece: C_H int_0^1 x^{N-5/2}(1+x/sqrt2)^{-N} dx <= sqrt(pi/2) 2^{3/4} * 2^{3/2} * exp(-0.5348*(N-5/2))
    T_0 = sqrt(pi/2)*2**(mpf(3)/4)*(2**(mpf(3)/2)*exp(-mpf('0.5348')*(N-mpf(5)/2)) + 2**(-N+mpf(5)/2)*2**(mpf(3)/2))
    q['Tails'] = T_R + T_L1 + T_0
    # ----- assemble
    bracket = E_mom + E_tailmom
    E_total = (2/pi)*q['eRe_phi']*q['H_up']*bracket + (2/pi)*q['Tails']
    q['E_total'] = E_total
    q['budget'] = P.C2*N**(-2)*exp(-th)
    return q

if __name__ == "__main__":
    mp.dps = 50
    for N in [400, 1000, 10**4, 10**5, 10**6, 10**8]:
        q = quant(N)
        ratio = q['E_total']/q['budget']
        print(f"N={N:9d}  E_total={float(q['E_total']):.3e}  budget={float(q['budget']):.3e}  ratio={float(ratio):.4f}")
        print(f"    E_mom={float(q['E_mom']):.3e}  E_tailmom={float(q['E_tailmom']):.3e}  Tails={float(q['Tails']):.3e}")
