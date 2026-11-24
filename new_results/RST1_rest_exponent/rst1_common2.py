"""RST1 v2 shared: generalized F_H lemma (rate lam), full assembled bound.
Lemma (Gaussian F_H, rate lam): for rL = sqrt(jq/(sqrt2 lam)) > c0 = 1+q/sqrt2:
F_H^{lam}(j,q) := int_0^inf e^{-lam r}|H| rho^j/|s+q|^{j+1} dr
 <= sqpi/q [ CL(lam) e^{-Sp/4} + ( (rL-c0)/(2 sq2) )^{-3/2} sqrt(2pi/c2L) + (sq2 rL)^{-3/2}/lam ]
     * e^{lam c0/2} e^{-Sp},   Sp = sqrt(2 sqrt2 lam j q),
CL(lam)=2^{3/2}/(3 lam), c2L = 8 sqrt2 j q/(4 rL+c0)^3.
Universal box (any j,q): F_H^{lam}(j,q) <= (2 sq2/q)(2 rL e^{lam c0/2}+1/lam) e^{-Sp}.
"""
import mpmath as mp
from mpmath import mpf, pi, sqrt, exp
mp.mp.dps = 60
j_=1j; omega=exp(j_*pi/4); a=pi/(2*sqrt(2)); sq2=sqrt(2); sqpi=sqrt(pi); sq2pi=sqrt(2*pi)
def H(s): return mp.beta(s/2, mpf(3)/2)/2

def FH_gauss(j,q,lam=a):
    rL=sqrt(j*q/(sq2*lam)); c0=1+q/sq2; Sp=sqrt(2*sq2*lam*j*q)
    assert rL>c0,(j,q,lam)
    CL=2**mpf('1.5')/(3*lam)
    c2=8*sq2*j*q/(4*rL+c0)**3
    mid=((rL-c0)/(2*sq2))**mpf('-1.5')*sqrt(2*pi/c2)
    return sqpi/q*(CL*exp(-Sp/4)+mid+(sq2*rL)**mpf('-1.5')/lam)*exp(lam*c0/2)*exp(-Sp)

def FH_box(j,q,lam=a):
    rL=sqrt(j*q/(sq2*lam)); c0=1+q/sq2; Sp=sqrt(2*sq2*lam*j*q)
    return 2*sq2/q*(2*rL*exp(lam*c0/2)+1/lam)*exp(-Sp)

def tg(k):
    N=k+1; return N**mpf('-0.5')*exp(-sq2pi*sqrt(N))

# --- pieces, each normalized to target tg(k) ---
def P_q2(k):
    return (8/pi)*(FH_gauss(k,2)+2*FH_gauss(k-1,2))/tg(k)
def P_mid(k):
    """i in [2,12], exact R2, each term normalized."""
    S=mpf(0)
    for i in range(2,13):
        S += (8/pi)*(2**i/mp.factorial(i))*FH_gauss(k-i,i+1)/tg(k)
    return S
def P_mid2(k):
    """i in [13,k/2]: box + convexity lump, exponent e^{-sqrt(14 pi (k-13))}."""
    if k < 26: return mpf(0)
    # lump: sum_i 2^i/i! * (2 sq2/(i+1))(2 rL e^{a c0/2}+1/a) with Q_i<=Q_13 (convexity)
    # prefactors bounded: rL <= sqrt(k(k+1)/pi)... use per-term to be safe but relax exponent:
    # exponent for i in [13,k/2]: (k-i)(i+1) >= 14(k-13) (min at i=13 for k/2>=13? min at endpoints)
    # we use the exact sum with box, exponent per term, but prefactor e^{a c0/2} bounded per term:
    S=mpf(0)
    for i in range(13, k//2+1):
        S += (8/pi)*(2**i/mp.factorial(i))*FH_box(k-i,i+1)/tg(k)
    return S
def P_tail(k):
    S=mpf(0)
    for i in range(k//2+1, k+1):
        S += (8/pi)*(2**i/mp.factorial(i))*FH_box(k-i,i+1)/tg(k)
    return S
def P_T2(k):
    return (4/pi)*FH_gauss(k,1,3*a)/tg(k)
def P_T3J(k):
    return (2/(pi*k)*3**(-k) + 2**(-k))/tg(k)

def totalC(k):
    return P_q2(k)+P_mid(k)+P_mid2(k)+P_tail(k)+P_T2(k)+P_T3J(k)

import re
def load_mk():
    out={}
    for line in open("<project>/research/08_theory/D13f_sweep/mk_deg9.tsv"):
        p=line.rstrip("\n").split("\t")
        if not p[0].isdigit(): continue
        out[int(p[0])]=mpf(re.match(r"\[([-0-9.eE+]+)",p[1]).group(1))
    return out
def Mmain(k):
    N=k+1; rs=mp.sqrt(2*N/pi)
    f=lambda r: exp(j_*pi*(r*omega)/2)*(r*omega/(1+r*omega))**N*H(r*omega)/(r*omega)*omega
    return 2/pi*mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]).imag
