import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")

# First: design the cell-based tail bounds using endpoint monotonicity.
# For a cell [Nlo, Nhi]:  y in [ylo, yhi] (y decreasing), kap^2 increasing.
# Tail pieces are written as products of factors, each bounded at a cell endpoint.
# Key quantities as functions of N (arb point evaluation needed only at endpoints).

newfile = '''"""dcy2_cell.py — cell-based certification via endpoint monotonicity.

For each piece of LHS = E_total * N^2 * exp(theta), every factor is explicit in
N (or y = sqrt(pi/(2N))), and on a cell [Nlo, Nhi] is bounded at an endpoint:
  y            : decreasing in N  -> use ylo for y^{+p}, yhi for y^{-p} or (1-y)^{-1}
  kap^2        : increasing       -> kap(Nhi) for kap^{+p}, kap(Nlo) for e^{-kap^2/4}
  kap*sqrt(y)  : decreasing for N >= 2956 -> max at Nlo (for N >= 2956)
  e^{-kap^2/4} = e^{-kap0^2/4} * (400/N) : decreasing -> at Nlo
  Delta = kap*sig : for N >= 2956, y*Delta = kap*sqrt(y)/sqrt(1.8*bR) decreasing
Pointwise (non-interval) arb evaluations at endpoints only; products certified.
"""
from flint import arb, acb, ctx

PI = arb.pi(); SQ2 = arb(2).sqrt(); bR = PI/(2*SQ2)
OMEGA = acb(1,1)/SQ2; WB = OMEGA.conjugate()
KAP0 = arb("2.48"); N1 = arb(400)
C_G = arb("1.0")

def y_of(N): return (PI/(2*arb(N))).sqrt()
def kap_of(N): return (KAP0**2 + 4*(arb(N)/N1).log()).sqrt()
def sig_of(N): return 1/(2*arb("0.9")*bR*y_of(N)).sqrt()
def Delta_of(N): return kap_of(N)*sig_of(N)

def ReL1max_of(N):
    y = y_of(N); return y*(1+bR*y)

def L_of(N, m):
    y = y_of(N); return arb(N)*y**(m+1) + y**m/m

def c2min_of(N): return arb("0.9")*bR*y_of(N)

def H_lo_of(N):
    y = y_of(N)
    return (PI/2).sqrt()*y**arb("1.5")*(1-3*y/4-C_G*y**2)
def H_up_of(N):
    y = y_of(N)
    return (PI/2).sqrt()*y**arb("1.5")*(1+3*y/4+C_G*y**2)

def alpha_of(N):
    D = Delta_of(N); L3,L4,L5 = L_of(N,3),L_of(N,4),L_of(N,5)
    return L3+L4*D+L5*D**2

def B3_of(N):
    y = y_of(N); D = Delta_of(N)
    return (arb(N)*y**4 + y**3/3)/(1-y*D)
def B6_of(N):
    y = y_of(N); D = Delta_of(N)
    return (arb(N)*y**7 + y**6/6)/(1-y*D)

def beta2_of(N):
    D = Delta_of(N); L3,L4,L5 = L_of(N,3),L_of(N,4),L_of(N,5)
    return L3*L5 + L4**2/2 + L4*L5*D + L5**2*D**2/2
def beta3_of(N):
    D = Delta_of(N); L3,L4,L5 = L_of(N,3),L_of(N,4),L_of(N,5)
    return (3*L3**2*L4 + (3*L3**2*L5+3*L3*L4**2)*D + (6*L3*L4*L5+L4**3)*D**2
            + (3*L3*L5**2+3*L4**2*L5)*D**3 + 3*L4*L5**2*D**4 + L5**3*D**5)/6

def M4p_of(N):
    y = y_of(N); D = Delta_of(N)
    x = (1/y-D)/SQ2
    return SQ2*(arb(11)/2).gamma()*x**(-arb(11)/2)/H_lo_of(N)

def binom(n, k):
    n = arb(n); k = int(k); r = arb(1)
    for i in range(1, k+1): r = r*(n-(i-1))/i
    return r

def Mtilde_of(N, jj, c2use=None):
    c2m = c2min_of(N) if c2use is None else c2use
    ReL1 = ReL1max_of(N); tR = ReL1/(2*c2m)
    tot = arb(0)
    for i in range(jj+1):
        tot += binom(jj,i)*tR**(jj-i)*c2m**(-(arb(i)+1)/2)*((arb(i)+1)/2).gamma()
    return (ReL1**2/(4*c2m)).exp()*tot

def eReL1Delta_of(N):
    return (ReL1max_of(N)*Delta_of(N)).exp()

def Gamma_upper_half(n, x):
    x = arb(x); res = PI.sqrt()*(x.sqrt()).erfc(); G = res
    for k in range(1, n+1):
        G = (k-arb(1)/2)*G + x**(k-arb(1)/2)*(-x).exp()
    return G
def Gamma_upper_int(n, x):
    x = arb(x); tot = arb(0)
    for i in range(n): tot += x**i/arb(i).fac()
    return arb(n-1).fac()*(-x).exp()*tot

def TailMom_of(N, p):
    c2m = c2min_of(N); s = (arb(p)+1)/2; z = c2m*Delta_of(N)**2/2
    G = Gamma_upper_int((p+1)//2, z) if p % 2 == 1 else Gamma_upper_half(p//2, z)
    return arb(1)/2*(c2m/2)**(-s)*G

def EdgeInt_of(N, jj, rate, T):
    s = (arb(jj)+1)/2; D = Delta_of(N)
    zT = rate*T**2; zD = rate*D**2
    if jj % 2 == 1:
        GT = Gamma_upper_int((jj+1)//2, zT); GD = Gamma_upper_int((jj+1)//2, zD)
    else:
        GT = Gamma_upper_half(jj//2, zT); GD = Gamma_upper_half(jj//2, zD)
    d = GT-GD
    if d < 0: d = arb(0)
    return arb(1)/2*rate**(-s)*d

def Gsplit_of(N, jj, gam):
    c2m = c2min_of(N); D = Delta_of(N)
    if gam*D <= c2m/2:
        return Mtilde_of(N, jj, c2m/2)
    T = min(D, c2m/(2*gam))
    base = (gam*T**3).exp()*Mtilde_of(N, jj)
    excess = gam*D - c2m
    e2 = eReL1Delta_of(N)
    if excess < 0:
        edge = EdgeInt_of(N, jj, -excess, T)*e2
    else:
        edge = e2*(excess*D**2).exp()*D**(jj+1)/(jj+1)
    return base + edge

def E_mom_of(N):
    q = {}
    y = y_of(N); D = Delta_of(N)
    Hsup = 1
    for jj in range(1,4):
        Hsup += [None, arb("1.5")*y*(1+2*y), arb("1.875")*y**2*(1+2*y), arb("2.1875")*y**3*(1+2*y)][jj]*D**jj
    L3,L4,L5 = L_of(N,3),L_of(N,4),L_of(N,5)
    B3,B6,al = B3_of(N), B6_of(N), alpha_of(N)
    be2, be3 = beta2_of(N), beta3_of(N)
    E = (M4p_of(N)/24)*Gsplit_of(N, 4, B3)
    E += Hsup*( be2*Mtilde_of(N,8) + be3*Mtilde_of(N,10) + al**4/24*Gsplit_of(N,12,al)
                + B6*(B6*D**6).exp()*Gsplit_of(N,6,al) )
    H1 = arb("1.5")*y*(1+2*y); H2 = arb("1.875")*y**2*(1+2*y); H3 = arb("2.1875")*y**3*(1+2*y)
    E_drop = {6: H1*L5+H2*L4+H3*L3, 7: H2*L5+H3*L4,
              8: H1*L3*L4+H2*L3**2/2+H3*L5, 9: H2*L3*L4+H3*L3**2/2,
              10: H1*L3**3/6, 11: H2*L3**3/6, 12: H3*L3**3/6}
    E += sum(cf*Mtilde_of(N,p) for p,cf in E_drop.items())
    return E

def E_tailmom_of(N):
    y = y_of(N); L3,L4,L5 = L_of(N,3),L_of(N,4),L_of(N,5)
    H1 = arb("1.5")*y*(1+2*y); H2 = arb("1.875")*y**2*(1+2*y); H3 = arb("2.1875")*y**3*(1+2*y)
    pp = {0: arb(1), 1: H1, 2: H2, 3: H3+L3, 4: L4+H1*L3, 5: L5+H1*L4+H2*L3,
          6: L3**2/2, 7: L3*L4+H1*L3**2/2, 9: L3**3/6}
    return sum(cf*TailMom_of(N,p) for p,cf in pp.items())

# ---- tail pieces via endpoint monotonicity (functions of N) ----
def S_sup(N):
    y = y_of(N); return ((PI*y/6)/(1-y)).exp()

def eLamD_up(N):
    """upper bound of exp(Lambda(Delta)) <= exp(-kap^2/4 + kap*c0*sqrt(y)), c0=0.79*sqrt2"""
    kap = kap_of(N); y = y_of(N)
    return (-kap**2/4 + kap*arb("1.117")*y.sqrt()).exp()

def CH_R_of(N):
    y = y_of(N); D = Delta_of(N)
    return (PI/2).sqrt()*2**(arb(3)/4)*y**arb("1.5")*(1+y*D)**(-arb("1.5"))

def rDelta_of(N):
    """lower bound of right-tail rate: bR*y*(2*Delta+y*Delta^2-0.514)/(1+y*Delta)^2"""
    y = y_of(N); D = Delta_of(N)
    return bR*y*(2*D+y*D**2-arb("0.514"))/(1+y*D)**2

def T_R_of(N):
    return (2/PI)*S_sup(N)*eLamD_up(N)*CH_R_of(N)/rDelta_of(N)*arb(N)**2

def LogM_of(N, t):
    y = y_of(N); u = 1+y*t
    Q = u**2 + SQ2*y*u + y**2
    return (arb(N)-1)*(u).log() + (arb(N)/2)*(1+SQ2*y+y**2).log() - (arb(N)/2)*(Q).log() - bR*t

def lamL_of(N):
    y = y_of(N); D = Delta_of(N); uD = 1-y*D
    QD = uD**2 + SQ2*y*uD + y**2
    N = arb(N)
    P = -2*bR*uD**3 - (PI+2)*y*uD**2 + (SQ2*(N-2)*y**2 - 2*bR*y**3)*uD + 2*(N-1)*y**3
    return min(P, PI*y/2)/(2*uD*QD)

def T_L1_of(N):
    y = y_of(N); D = Delta_of(N)
    CH_L = (PI/2).sqrt()*2**(arb(3)/4)
    return (2/PI)*S_sup(N)*LogM_of(N,-D).exp()*CH_L/lamL_of(N)*arb(N)**2

def T_0_of(N):
    N = arb(N)
    th = (PI*arb(N)).sqrt()
    return (2/PI)*(PI/2).sqrt()*2**(arb(3)/4)*(2**(arb(3)/2)*(-arb("0.5348")*(N-arb(5)/2)).exp()
        + 2**(-N+arb(5)/2)*2**(arb(3)/2))*N**2*th.exp()

def LHS_of(N):
    """E_total*N^2*exp(theta) = (2/pi)*y*e^{S_sup}*[H_up*bracket + tailpieces]*N^2 + T_0"""
    y = y_of(N)
    bracket = E_mom_of(N) + E_tailmom_of(N)
    tailpieces = eLamD_up(N)*CH_R_of(N)/rDelta_of(N) + LogM_of(N,-Delta_of(N)).exp()*(PI/2).sqrt()*2**(arb(3)/4)/lamL_of(N)
    return (2/PI)*y*S_sup(N)*(H_up_of(N)*bracket + tailpieces)*arb(N)**2 + T_0_of(N)
'''
open("dcy2_cell.py","w").write(newfile)
import subprocess
test = '''
from flint import ctx
ctx.dps = 40
import dcy2_cell as C
C2 = 10**6
for N in [400, 2000, 3000, 10000, 10**5, 10**6, 10**9, 10**12]:
    v = C.LHS_of(N)
    print(f"N={N:9d}  LHS={v.str(10)}   LHS/C2={float(v/C2):.5f}")
'''
open("t_lhs.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_lhs.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
