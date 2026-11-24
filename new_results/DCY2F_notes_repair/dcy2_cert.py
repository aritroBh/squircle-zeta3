"""dcy2_cert.py — arb-certified bounds for DCY2 (outward rounding throughout).

All quantities are functions of N (arb, possibly an interval [Nlo,Nhi]).
The structure mirrors dcy2_budget.py but every bound is a certified arb ball.

Notation: N >= 400, y = sqrt(pi/(2N)), rho = 1/y, theta = sqrt(pi N) = pi/(sqrt2 y),
omega = e^{i pi/4}, wb = conjugate, b = pi wb/2, bR = Re b = pi/(2 sqrt2).
"""
from flint import arb, acb, ctx

# ---------------- helpers ----------------
def iv(lo, hi):
    """proper endpoint interval [lo,hi] with outward rounding."""
    lo = arb(lo); hi = arb(hi)
    mid = (lo+hi)/2
    rad = (hi-lo)/2
    return arb(mid, rad)

def I(x):
    return arb(x)

def binom(n, k):
    n = arb(n); k = int(k)
    r = arb(1)
    for i in range(1, k+1):
        r = r*(n-(i-1))/i
    return r

def IC(x):
    return acb(x)

PI = arb.pi()
SQ2 = arb(2).sqrt()
OMEGA = acb(1, 1)/SQ2      # e^{i pi/4}
WB = OMEGA.conjugate()
bR = PI/(2*SQ2)
bC = PI*WB/2               # complex b (for reference)
KAP0 = arb("2.48")
N1 = arb(400)

# ---------------- Watson / H constants (certified in step_03) ----------------
C_G   = arb("1.0")      # |g(s0)-1+3/(4 s0)| <= C_G y^2  (certified in step_05_watson_g)
C_H1R = arb("1.0")      # |h1| <= (3/2) y (1+2y)  [Watson form]
C_H2R = arb("1.0")      # |h2 - (15/8) y^2| <= C_H2R y^3
C_H3R = arb("2.5")      # |h3| <= C_H3R y^3  (= (35/16)(1+2 y1) + rem)

# ---------------- core quantities ----------------
def quantities(N):
    q = {}
    N = arb(N)
    q['N'] = N
    q['y'] = y = (PI/(2*N)).sqrt()
    q['rho'] = rho = 1/y
    q['theta'] = theta = PI/(SQ2*y)
    q['kap2'] = kap2 = KAP0**2 + 4*(N/N1).log()
    q['kap'] = kap = kap2.sqrt()
    q['s0'] = s0 = OMEGA/y
    # --- c2 = -Re l2 >= 0.9 bR y (Lemma E2; certified poly inequality in step_04)
    q['c2min'] = c2m = arb("0.9")*bR*y
    q['sig'] = sig = 1/(2*c2m).sqrt()
    q['Delta'] = Delta = kap*sig
    q['yDelta'] = yDelta = y*Delta
    # --- |Re l1| <= y(1+bR y) (exact formula, Lemma E1)
    q['ReL1max'] = y*(1+bR*y)
    # --- |l_m| <= L_m for m >= 2 (Lemma E3)
    def L(m):
        return N*y**(m+1) + y**m/m
    q['L'] = L
    # --- e^{Re phi(s0)} <= y e^{-theta} e^{S_sup}, S_sup = (pi y/6)/(1-y)
    S_sup = (PI*y/6)/(1-y)
    q['eRe_phi'] = y*(-theta).exp()*S_sup.exp()
    # --- |H(s0)| bounds: |g| in [1-3y/4-C_G y^2, 1+3y/4+C_G y^2]
    g_lo = 1 - 3*y/4 - C_G*y**2
    g_up = 1 + 3*y/4 + C_G*y**2
    q['H_lo'] = (PI/2).sqrt()*y**arb("1.5")*g_lo
    q['H_up'] = (PI/2).sqrt()*y**arb("1.5")*g_up
    # --- |h_j| bounds: H_j = (3/2)_j y^j/j! (1+2y)  [valid when C_HjR y <= 2 - j/2]
    def dblfac_ratio(jj):
        # (3/2)_j / j! = (2j+1)!!/(2^j j!)
        r = arb(1)
        for k in range(1, jj+1):
            r *= (arb(k)+arb(1)/2)/k
        return r
    q['H'] = [arb(1)] + [dblfac_ratio(jj)*y**jj*(1+2*y) for jj in range(1, 4)]
    # --- series sups
    B3 = (N*y**4 + y**3/3)/(1-yDelta)
    B6 = (N*y**7 + y**6/6)/(1-yDelta)
    q['B3'] = B3; q['B6'] = B6
    L3, L4, L5 = L(3), L(4), L(5)
    alpha = L3 + L4*Delta + L5*Delta**2
    q['alpha'] = alpha
    beta2 = L3*L5 + L4**2/2 + L4*L5*Delta + L5**2*Delta**2/2
    beta3 = (3*L3**2*L4 + (3*L3**2*L5+3*L3*L4**2)*Delta + (6*L3*L4*L5+L4**3)*Delta**2
             + (3*L3*L5**2+3*L4**2*L5)*Delta**3 + 3*L4*L5**2*Delta**4 + L5**3*Delta**5)/6
    # V-factors
    q['V3'] = V3 = max(arb(0), B3*Delta**3 - c2m*Delta**2)
    q['VA'] = VA = max(arb(0), alpha*Delta**3 - c2m*Delta**2)
    # --- M4' (H fourth-moment remainder, valid for all t in window)
    x = (1/y-Delta)/SQ2
    M4num = SQ2*(arb(11)/2).gamma()*x**(-arb(11)/2)
    q['M4'] = M4num / q['H_lo']
    # --- |cdot| moments (Gaussian kept): Mtilde_j
    ReL1 = q['ReL1max']
    tR = ReL1/(2*c2m)
    efac = (ReL1**2/(4*c2m)).exp()
    def Mtilde(jj, c2use=None):
        cc = c2m if c2use is None else c2use
        tR2 = ReL1/(2*cc)
        tot = arb(0)
        for i in range(jj+1):
            tot += binom(jj, i) * tR2**(jj-i) * cc**(-(arb(i)+1)/2) * ((arb(i)+1)/2).gamma()
        return (ReL1**2/(4*cc)).exp()*tot
    q['Mtilde'] = Mtilde
    # --- plain moments (no Gaussian): int_{|t|<=Delta} |t|^j e^{|Re l1| |t|} dt
    def PMtilde(jj):
        return 2*(ReL1*Delta).exp()*Delta**(jj+1)/(jj+1)
    q['PMtilde'] = PMtilde
    def Gamma_upper_half(n, x):
        # Gamma(n+1/2, x) for integer n >= 0, via recurrence from Gamma(1/2,x) = sqrt(pi) erfc(sqrt x)
        res = PI.sqrt()*(x.sqrt()).erfc()
        out = {0: res}
        G = res
        for k in range(1, n+1):
            G = (k-arb(1)/2)*G + x**(k-arb(1)/2)*(-x).exp()
            out[k] = G
        return out[n]
    def Gamma_upper_int(n, x):
        # Gamma(n, x) = (n-1)! e^{-x} sum_{i=0}^{n-1} x^i/i!
        tot = arb(0)
        for i in range(n):
            tot += x**i/arb(i).fac()
        return arb(n-1).fac()*(-x).exp()*tot

    # --- moment with e^{gamma |t|^3} factor:
    # base: e^{gamma T^3} * Mtilde (full line), T = min(Delta, c2m/(2 gamma));
    # edge [T,Delta]: int |t|^j e^{(gamma Delta - c2m) t^2} dt, damped Gaussian tail.
    def EdgeInt(jj, rate, T):
        # int_T^Delta |t|^j e^{-rate * t^2} dt  (rate > 0):
        # = (1/2) rate^{-(j+1)/2} [Gamma((j+1)/2, rate T^2) - Gamma((j+1)/2, rate Delta^2)]
        s = (arb(jj)+1)/2
        zT = rate*T**2
        zD = rate*Delta**2
        if jj % 2 == 1:
            GT = Gamma_upper_int((jj+1)//2, zT)
            GD = Gamma_upper_int((jj+1)//2, zD)
        else:
            GT = Gamma_upper_half(jj//2, zT)
            GD = Gamma_upper_half(jj//2, zD)
        d = GT - GD
        if d < 0:
            d = arb(0)
        return arb(1)/2*rate**(-s)*d
    def Gsplit(jj, gam):
        if gam*Delta <= c2m/2:
            # gamma |t|^3 <= (c2m/2) t^2 on the whole window: widened Gaussian moment
            return Mtilde(jj, c2m/2)
        T = min(Delta, c2m/(2*gam))
        base = (gam*T**3).exp()*Mtilde(jj)
        excess = gam*Delta - c2m
        efac2 = (ReL1*Delta).exp()
        if excess < 0:
            edge = EdgeInt(jj, -excess, T)*efac2
        else:
            edge = efac2*(excess*Delta**2).exp()*Delta**(jj+1)/(jj+1)
        return base + edge
    q['Gsplit'] = Gsplit
    # --- H sup on window: |1+h1 t+h2 t^2+h3 t^3| <= Hsup
    Hsup = 1 + q['H'][1]*Delta + q['H'][2]*Delta**2 + q['H'][3]*Delta**3
    q['Hsup'] = Hsup
    # --- E_mom
    eV3 = V3.exp(); eVA = VA.exp(); eB6 = (B6*Delta**6).exp()
    E_mom = (q['M4']/24)*Gsplit(4, B3)
    E_mom += Hsup*( beta2*Mtilde(8) + beta3*Mtilde(10) + alpha**4/24*Gsplit(12, alpha)
                    + B6*eB6*Gsplit(6, alpha) )
    E_drop = {
        6:  q['H'][1]*L5 + q['H'][2]*L4 + q['H'][3]*L3,
        7:  q['H'][2]*L5 + q['H'][3]*L4,
        8:  q['H'][1]*L3*L4 + q['H'][2]*L3**2/2 + q['H'][3]*L5,
        9:  q['H'][2]*L3*L4 + q['H'][3]*L3**2/2,
        10: q['H'][1]*L3**3/6,
        11: q['H'][2]*L3**3/6,
        12: q['H'][3]*L3**3/6,
    }
    E_mom += sum(cf*Mtilde(p) for p, cf in E_drop.items())
    q['E_mom'] = E_mom
    # --- window-tail moments (Gaussian tail beyond Delta)
    def TailMom(p):
        s = (arb(p)+1)/2
        z = c2m*Delta**2/2
        if p % 2 == 1:
            G = Gamma_upper_int((p+1)//2, z)
        else:
            G = Gamma_upper_half(p//2, z)
        return arb(1)/2*(c2m/2)**(-s)*G
    q['TailMom'] = TailMom
    pp = {
        0: arb(1),
        1: q['H'][1],
        2: q['H'][2],
        3: q['H'][3] + L3,
        4: L4 + q['H'][1]*L3,
        5: L5 + q['H'][1]*L4 + q['H'][2]*L3,
        6: L3**2/2,
        7: L3*L4 + q['H'][1]*L3**2/2,
        9: L3**3/6,
    }
    q['E_tailmom'] = sum(cf*TailMom(p) for p, cf in pp.items())
    # --- outer tails
    def LogM(t):
        u = 1+y*t
        Q = u**2 + SQ2*y*u + y**2
        return (N-1)*(u).log() + (N/2)*(1+SQ2*y+y**2).log() - (N/2)*(Q).log() - bR*t
    LM_D = LogM(Delta)
    LM_mD = LogM(-Delta)
    # right tail rate
    uDR = 1+yDelta
    QR = uDR**2 + SQ2*y*uDR + y**2
    LamPR = y*((N-1)/uDR - N*(uDR+y/SQ2)/QR) - bR
    rDelta = min(-LamPR, bR)
    lamR = rDelta
    # left rate
    uD = 1-yDelta
    QD = uD**2 + SQ2*y*uD + y**2
    def Pofu(u):
        return -2*bR*u**3 - (PI+2)*y*u**2 + (SQ2*(N-2)*y**2 - 2*bR*y**2)*u + 2*(N-1)*y**3
    lamL = min(Pofu(uD), PI*y/2)/(2*uD*QD)
    q['lamL'] = lamL
    CH_R = (PI/2).sqrt()*((rho+Delta)/SQ2)**(-arb(3)/2)
    CH_L = (PI/2).sqrt()*2**(arb(3)/4)
    T_R = q['eRe_phi']*LM_D.exp()*CH_R/lamR
    T_L1 = q['eRe_phi']*LM_mD.exp()*CH_L/lamL
    T_0 = (PI/2).sqrt()*2**(arb(3)/4)*(2**(arb(3)/2)*(-arb("0.5348")*(N-arb(5)/2)).exp()
          + 2**(-N+arb(5)/2)*2**(arb(3)/2))
    q['Tails'] = T_R + T_L1 + T_0
    # --- assemble
    bracket = E_mom + q['E_tailmom']
    q['E_total'] = (2/PI)*q['eRe_phi']*q['H_up']*bracket + (2/PI)*q['Tails']
    # --- LHS = E_total * N^2 * e^{theta}, with the e^{-theta} factors cancelled algebraically.
    # eRe_phi = y e^{-theta} e^{S_sup} is folded into [H_up*bracket + tail-pieces].
    S_sup = (PI*y/6)/(1-y)
    tailpieces = LM_D.exp()*CH_R/lamR + LM_mD.exp()*CH_L/lamL
    LHS = (2/PI)*y*S_sup.exp()*( q['H_up']*bracket + tailpieces )*N**2           + (2/PI)*T_0*N**2*q['theta'].exp()
    q['LHS'] = LHS
    return q

C2 = arb(10)**6

def budget(N):
    q = quantities(N)
    return C2*q['N']**(-2)*(-q['theta']).exp()

def ratio(N):
    q = quantities(N)
    return q['LHS']/C2, q

if __name__ == "__main__":
    ctx.dps = 40
    for N in [400, 2000, 10000, 10**6]:
        r, q = ratio(N)
        print(f"N={N:10d}  ratio={r.str(12)}")
