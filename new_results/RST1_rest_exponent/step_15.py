"""RST1 step_15: arb-certified constants for the theorem.
Every piece is non-increasing for k>=99 (monotonicity, hand proof in proof.md);
the constant is the sum of the arb-evaluated pieces at k=99.
"""
from flint import arb, acb, ctx
ctx.dps = 80

def A(x): return arb(x)
pi = arb.pi(); sq2 = arb(2).sqrt(); sqpi = pi.sqrt(); sq2pi = (2*pi).sqrt()
a = pi/(2*sq2)
def fac(n): return arb(n).fac()

# ---- F_H bounds (arb) ----
def FH_gauss(j, q, lam):
    j, q = arb(j), arb(q)
    rL = (j*q/(sq2*lam)).sqrt(); c0 = 1 + q/sq2
    Sp = (2*sq2*lam*j*q).sqrt()
    assert rL > c0
    CL = 2**arb('1.5')/(3*lam)
    c2 = 8*sq2*j*q/(4*rL+c0)**3
    mid = ((rL-c0)/(2*sq2))**arb('-1.5') * (2*pi/c2).sqrt()
    return sqpi/q*(CL*(-Sp/4).exp() + mid + (sq2*rL)**arb('-1.5')/lam)*(lam*c0/2).exp()*(-Sp).exp()

def FH_box(j, q, lam):
    j, q = arb(j), arb(q)
    rL = (j*q/(sq2*lam)).sqrt(); c0 = 1 + q/sq2
    Sp = (2*sq2*lam*j*q).sqrt()
    return 2*sq2/q*(2*rL*(lam*c0/2).exp() + 1/lam)*(-Sp).exp()

def tg(k):
    N = arb(k+1)
    return N**arb('-0.5')*(-sq2pi*N.sqrt()).exp()

K = 99
# q2 (i=0, i=1)
q2 = 8/pi*(FH_gauss(K,2,a) + 2*FH_gauss(K-1,2,a))/tg(K)
# mid (i in [2,12])
mid = sum(8/pi*(arb(2)**i/fac(i))*FH_gauss(K-i,i+1,a)/tg(K) for i in range(2,13))
# mid2 (i in [13,K/2]): box + convexity(Q_i <= Q_13) lump
if K >= 26:
    xsum = sum(arb(2)**i/fac(i) for i in range(13, K//2+1))
    c013 = 1 + 14/sq2
    Q13 = (a*c013/2).exp() * (-(14*pi*(K-13)).sqrt()).exp()
    rLm = (K*(K+1)/pi).sqrt()
    mid2 = 8/pi * 2*sq2/14 * (2*rLm*Q13 + (1/a)*(-(14*pi*(K-13)).sqrt()).exp()) * xsum / tg(K)
else:
    mid2 = arb(0)
# tail (i > K/2): drop exponent, CX6(16) coefficient trick
xp = 2*(a/(2*sq2)).exp()
tail_sum1 = (a/(2*sq2)).exp() * (arb(1).exp()*xp).exp()   # e^{a/(2 sqrt2)} e^{e x'}
tail_sum2 = (2*arb(1).exp()).exp()                          # e^{2e}
rLt = (K*(K+1)/pi).sqrt()
tail = 8/pi * 2*sq2/(K//2+1) * (2*rLt*tail_sum1 + tail_sum2/a) * (-arb(K)/2).exp() / tg(K)
# T2
T2 = 4/pi * FH_gauss(K,1,3*a)/tg(K)
# T3J
T3J = (2/(pi*K)*arb(3)**(-K) + arb(2)**(-K))/tg(K)

C_raw = q2 + mid + mid2 + tail + T2 + T3J
for nm, v in [('a',a),('sqrt(2 pi)',sq2pi),('q2',q2),('mid',mid),('mid2',mid2),('tail',tail),('T2',T2),('T3J',T3J),('C_raw',C_raw)]:
    print(f"{nm} = {v.str(22)}", flush=True)
C = 1030
print(f"C = {C}; margin C - C_raw = {(C - C_raw).str(10)}", flush=True)
assert C_raw < C and (C - C_raw)/C_raw > 0.001
print("CERTIFIED: C_raw < 1030", flush=True)
