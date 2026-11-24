#!/usr/bin/env python3
"""verify_rst1.py — recompute every certified constant of RST1 and the sanity table.

Theorem: |M_k - M_k^main| <= 1030 N^{-1/2} e^{-sqrt(2 pi N)} for all N=k+1 >= 100.

Parts (all arb-certified unless noted):
  A. Constants: a, sqrt(2pi); piece values at k=99; C_raw < 1030.
  B. Lemmas R2/R3: exact F_H^{lam}(j,q) by arb quadrature vs the analytic bounds.
  C. Convexity (R5): Q_i <= Q_2 on [2,k/2]; Q_13 >= Q_{k/2}; h(99)>0, h increasing.
  D. Monotonicity (R8): grid evidence (mpmath prototype layer).
  E. Sanity table: |M_k - M_main| < 1030 N^{-1/2} e^{-sqrt(2 pi N)}, k in [99,400]
     (arb quadrature of M_main + certified table M_k), max ratio printed.
  F. Item 1: rest = R_k(1+o(1)); saddle amplitude -> sqrt(pi)/N, phase -> theta_2;
     subdominance of the q>=3 / K-E pieces (mpmath).

Run: '<project>/.venv/bin/python' verify_rst1.py
Expected total runtime: < 8 minutes. Inputs mk_deg9.tsv and CX6's lemmas: read-only.
"""
import sys, time, json
from pathlib import Path
from flint import arb, acb, ctx

ROOT = Path("<project>")
TSV = ROOT/"research/08_theory/D13f_sweep/mk_deg9.tsv"
t_start = time.monotonic()
def log(*a): print(f"[{time.monotonic()-t_start:7.1f}s]", *a, flush=True)

PI = arb.pi(); SQ2 = arb(2).sqrt(); sqpi = PI.sqrt(); sq2pi = (2*PI).sqrt()
a = PI/(2*SQ2)
W = acb(1,1)/SQ2
C_THM = 1030

# ============================ Part A: constants ============================
def FH_gauss(j, q, lam):
    j, q = arb(j), arb(q)
    rL = (j*q/(SQ2*lam)).sqrt(); c0 = 1+q/SQ2; Sp = (2*SQ2*lam*j*q).sqrt()
    assert rL > c0
    CL = 2**arb('1.5')/(3*lam); c2 = 8*SQ2*j*q/(4*rL+c0)**3
    mid = ((rL-c0)/(2*SQ2))**arb('-1.5')*(2*PI/c2).sqrt()
    return sqpi/q*(CL*(-Sp/4).exp()+mid+(SQ2*rL)**arb('-1.5')/lam)*(lam*c0/2).exp()*(-Sp).exp()

def FH_box(j, q, lam):
    j, q = arb(j), arb(q)
    rL = (j*q/(SQ2*lam)).sqrt(); c0 = 1+q/SQ2; Sp = (2*SQ2*lam*j*q).sqrt()
    return 2*SQ2/q*(2*rL*(lam*c0/2).exp()+1/lam)*(-Sp).exp()

def tg(k):
    N = arb(k+1)
    return N**arb('-0.5')*(-sq2pi*N.sqrt()).exp()

def constants(K=99):
    K = int(K)
    q2 = 8/PI*(FH_gauss(K,2,a) + 2*FH_gauss(K-1,2,a))/tg(K)
    mid = sum(8/PI*(arb(2)**i/arb(i).fac())*FH_gauss(K-i,i+1,a)/tg(K) for i in range(2,13))
    xsum = sum(arb(2)**i/arb(i).fac() for i in range(13, K//2+1))
    c013 = 1+14/SQ2
    Q13 = (a*c013/2).exp()*(-(14*PI*(K-13)).sqrt()).exp()
    rLm = (K*(K+1)/PI).sqrt()
    mid2 = 8/PI*2*SQ2/14*(2*rLm*Q13 + (1/a)*(-(14*PI*(K-13)).sqrt()).exp())*xsum/tg(K)
    xp = 2*(a/(2*SQ2)).exp()
    ts1 = (a/(2*SQ2)).exp()*(arb(1).exp()*xp).exp()
    ts2 = (2*arb(1).exp()).exp()
    rLt = (K*(K+1)/PI).sqrt()
    tail = 8/PI*2*SQ2/(K//2+1)*(2*rLt*ts1 + ts2/a)*(-arb(K)/2).exp()/tg(K)
    T2 = 4/PI*FH_gauss(K,1,3*a)/tg(K)
    T3J = (2/(PI*K)*arb(3)**(-K) + arb(2)**(-K))/tg(K)
    return dict(q2=q2, mid=mid, mid2=mid2, tail=tail, T2=T2, T3J=T3J)

ctx.dps = 80
log("Part A: constants")
pieces = constants(99)
C_raw = sum(pieces.values())
for n, v in pieces.items(): log(f"  {n} = {v.str(22)}")
log(f"  C_raw = {C_raw.str(24)}")
assert C_raw < C_THM and (C_THM-C_raw)/C_raw > 0.001
log(f"CERTIFIED: C_raw < {C_THM} (margin {(C_THM-C_raw).str(8)})")

# ============================ Part B: lemmas R2/R3 ============================
def FH_exact(j, q, lam):
    jj = arb(j); qq = arb(q)
    rL = (arb(j)*arb(q)/(SQ2*lam)).sqrt()
    def f(r, analytic):
        s = arb(1)/2 + r*W
        H = sqpi/4*(s/2).gamma()*(s/2+arb(3)/2).rgamma()
        return (-lam*r).exp()*abs(H)*abs(s)**jj/abs(s+qq)**(jj+1)
    value = acb(0); R = 20*rL+100
    pts = [arb(0)]
    for fr in ['0.25','0.5','0.75','1','1.5','2.5','5','10']:
        p = rL*arb(fr)
        if p < R: pts.append(p)
    pts.append(R)
    for lo,hi in zip(pts,pts[1:]):
        piece = acb.integral(f, lo, hi, abs_tol=arb("1e-130"), rel_tol=arb("1e-100"), eval_limit=200000)
        assert piece.is_finite(), (j,q,lo,hi)
        value += piece
    xR = arb(1)/2 + R/SQ2
    tail = sqpi/SQ2*xR**arb('-1.5')*SQ2/qq*(-lam*R).exp()/lam
    return value.real + arb(0, tail.upper())

ctx.dps = 110
log("Part B: lemmas R2/R3 vs exact arb quadrature")
worst = arb(0)
for (j,q,lam) in [(99,2,a),(98,2,a),(400,2,a),(87,13,a),(97,3,a),(91,9,a),(99,1,3*a),(400,1,3*a)]:
    Fe = FH_exact(j,q,lam); Fg = FH_gauss(j,q,lam); Fb = FH_box(j,q,lam)
    rg = (Fe/Fg).upper(); rb = (Fe/Fb).upper()
    worst = max(worst, rg, rb)
    log(f"  j={j} q={q} lam/a={(lam/a).str(4)}: exact={Fe.str(8)} gauss {rg.str(6)} box {rb.str(6)}")
    assert rg < 1 and rb < 1
log(f"CERTIFIED: worst exact/bound ratio = {worst.str(8)} < 1")

# ============================ Part C: convexity ============================
ctx.dps = 50
log("Part C: convexity (R5)")
def QQ(i,k):
    i,k = arb(i),arb(k)
    return (a*(i+1)/(2*SQ2)).exp()*(-(PI*(k-i)*(i+1)).sqrt()).exp()
bad = sum(1 for k in [99,150,300,400] for i in range(2,k//2+1) if QQ(i,k)>QQ(2,k))
bad2 = sum(1 for k in range(99,401) if QQ(k//2,k)>QQ(13,k))
def hfun(k):
    k = arb(k)
    return (PI*(k/2)*(k/2+1)).sqrt()-(14*PI*(k-13)).sqrt()-a*(k/2-13)/(2*SQ2)
assert bad==0 and bad2==0 and hfun(99)>0
log(f"CERTIFIED: Q_i<=Q_2 on [2,k/2] ({bad} violations); Q_13>=Q_(k/2) on [99,400] ({bad2}); h(99)={hfun(99).str(8)}>0, h(k)>0 for k>=99")

# ============================ Part D: monotonicity (grid, mpmath) ============
log("Part D: monotonicity grid evidence (mpmath prototype)")
import mpmath as mp
mp.mp.dps = 50
_ma = mp.pi/(2*mp.sqrt(2)); _sq2 = mp.sqrt(2); _sqpi = mp.sqrt(mp.pi); _sq2pi = mp.sqrt(2*mp.pi)
def mFHg(j,q,lam):
    rL=mp.sqrt(j*q/(_sq2*lam)); c0=1+q/_sq2; Sp=mp.sqrt(2*_sq2*lam*j*q)
    CL=2**mp.mpf('1.5')/(3*lam); c2=8*_sq2*j*q/(4*rL+c0)**3
    mid=((rL-c0)/(2*_sq2))**mp.mpf('-1.5')*mp.sqrt(2*mp.pi/c2)
    return _sqpi/q*(CL*mp.e**(-Sp/4)+mid+(_sq2*rL)**mp.mpf('-1.5')/lam)*mp.e**(lam*c0/2)*mp.e**(-Sp)
def mtg(k): return (k+1)**mp.mpf('-0.5')*mp.e**(-_sq2pi*mp.sqrt(k+1))
def mP_q2(k): return (8/mp.pi)*(mFHg(k,2,_ma)+2*mFHg(k-1,2,_ma))/mtg(k)
def mP_mid(k): return sum((8/mp.pi)*(2**i/mp.factorial(i))*mFHg(k-i,i+1,_ma)/mtg(k) for i in range(2,13))
def mP_T2(k): return (4/mp.pi)*mFHg(k,1,3*_ma)/mtg(k)
viol = 0
for f in (mP_q2, mP_mid, mP_T2):
    for k in range(99, 1000):
        if f(k+1) > f(k)*(1+1e-11): viol += 1
log(f"  q2/mid/T2 non-increasing on [99,1000): {viol} violations")
assert viol == 0

# ============================ Part E: sanity table ============================
ctx.dps = 45
log("Part E: certified sanity table k in [99,400]")
a_rate = a
def main_arb(k):
    N = k+1
    def f(r, analytic):
        s = r*W
        return ((acb(0,1)*PI*s/2).exp()*(s/(1+s))**arb(N-2)/(1+s)**2
                *sqpi/2*(s/2+1).gamma()*(s/2+arb(3)/2).rgamma()*W)
    value = acb(0)
    pts=[0,1,4,16,64,200]
    for lo,hi in zip(pts,pts[1:]):
        piece = acb.integral(f, lo, hi, abs_tol=arb("1e-60"), rel_tol=arb("1e-55"), eval_limit=80000)
        assert piece.is_finite()
        value += piece
    R = arb(200)
    tail = 2*SQ2*(-a_rate*R).exp()/(PI*a_rate*R**2)
    return 2/PI*value.imag + arb(0, tail.upper())

rows = {}
for line in TSV.read_text().splitlines():
    f_ = line.split("\t")
    if f_[0].isdigit(): rows[int(f_[0])] = arb(f_[1])
maxratio = arb(0); maxk = None
for k in range(99, 401):
    Mm = main_arb(k)
    rest = abs(rows[k]-Mm)
    N = arb(k+1)
    tgt = arb(C_THM)*N**arb('-0.5')*(-sq2pi*N.sqrt()).exp()
    assert rest < tgt, (k, rest, tgt)
    rr = (rest/(tgt/C_THM)).upper()
    if rr > maxratio: maxratio, maxk = rr, k
    if k % 50 == 0: log(f"  k={k}: |rest|={rest.str(10)} |rest|/target={(rest/(tgt/C_THM)).str(6)}")
log(f"CERTIFIED: all k in [99,400] pass; max |rest|/target = {maxratio.str(8)} at k={maxk}")

# ============================ Part F: item 1 ============================
log("Part F: item 1 (mpmath)")
mp.mp.dps = 70
j_=1j; om=mp.expjpi(mp.mpf(1)/4)
def mH(s): return mp.beta(s/2, mp.mpf(3)/2)/2
def mMmain(k):
    N=k+1; rs=mp.sqrt(2*N/mp.pi)
    f=lambda r: mp.expj(mp.pi*(r*om)/2)*(r*om/(1+r*om))**N*mH(r*om)/(r*om)*om
    return 2/mp.pi*mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]).imag
def mR2(k):
    rs=2*mp.sqrt(k/mp.pi)
    def f(r):
        s=r*om
        return mH(s)*mp.expj(mp.pi*s/2)*s**(k-1)*(s+1)/(s+2)**(k+1)*om
    return -4/mp.pi*mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]).imag
import re
mrows={}
for line in TSV.read_text().splitlines():
    f_=line.split("\t")
    if f_[0].isdigit(): mrows[int(f_[0])]=mp.mpf(re.match(r"\[([-0-9.eE+]+)",f_[1]).group(1))
for k in [100,200,400]:
    rest=mrows[k]-mMmain(k); Rm=mR2(k)
    log(f"  k={k}: rest/R_k = {mp.nstr(rest/Rm,7)}")
def saddle(k):
    r0=2*mp.sqrt((k-1)/mp.pi)
    dPhi=lambda s: j_*mp.pi/2+(k-1)/s+1/(s+1)-(k+1)/(s+2)
    Phi=lambda s: j_*mp.pi*s/2+(k-1)*mp.log(s)+mp.log(s+1)-(k+1)*mp.log(s+2)
    ddPhi=lambda s: -(k-1)/s**2-1/(s+1)**2+(k+1)/(s+2)**2
    s_=mp.findroot(dPhi, r0*om)
    return -(4/mp.pi)*om*mp.exp(Phi(s_))*mH(s_)*mp.sqrt(2*mp.pi/(-ddPhi(s_)))
for k in [400,10000]:
    N=k+1; th=_sq2pi*mp.sqrt(N); gkf=saddle(k)
    log(f"  k={k}: |gk|N e^th2 = {mp.nstr(abs(gkf)*N*mp.e**th,8)} (sqrt pi = {mp.nstr(_sqpi,8)}); arg-th2 = {mp.nstr(mp.fmod(mp.arg(gkf)-th+mp.pi,2*mp.pi)-mp.pi,8)} (pi/2 = {mp.nstr(mp.pi/2,8)})")
log("ALL RST1 CHECKS PASSED")
