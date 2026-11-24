"""RST1R step_04: R4 T1 assembly audit.
(a) coverage: q2(i=0,1) + mid[2,12] + mid2[13,k//2] + tail[k//2+1,k] covers i=0..k of (12) exactly once.
(b) mid2 (13) vs the exact box-sum at k=99 and other k (arb).
(c) tail (14) vs the exact box-sum at k=99 (arb) -- test the suspected missing e^{a/2}.
(d) per-term prefactor bounds used in (13)/(14): 2sqrt2/(i+1) <= 2sqrt2/14 resp <= 2sqrt2/(k/2+1);
    r_a <= sqrt(k(k+1)/pi); (k-i)(i+1) >= 14(k-13) for integer i in [13,k//2], all k>=26.
"""
from flint import arb, ctx
ctx.dps = 80
PI = arb.pi(); SQ2 = arb(2).sqrt(); sqpi = PI.sqrt()
a = PI/(2*SQ2)

def FH_box(j,q,lam=a):
    j,q = arb(j),arb(q)
    rL=(j*q/(SQ2*lam)).sqrt(); c0=1+q/SQ2; Sp=(2*SQ2*lam*j*q).sqrt()
    return 2*SQ2/q*(2*rL*(lam*c0/2).exp()+1/lam)*(-Sp).exp()

print("== (a) coverage ==")
for k in [99, 100, 101, 400]:
    idx = [0,1] + list(range(2,13)) + list(range(13,k//2+1)) + list(range(k//2+1,k+1))
    assert idx == list(range(0,k+1)), k
print("  split covers i=0..k exactly once: OK (k=99,100,101,400)")

print("== (b) mid2: exact box-sum vs (13) ==")
def mid2_exact(k):
    return sum(arb(2)**i/arb(i).fac()*FH_box(k-i,i+1) for i in range(13, k//2+1))
def mid2_bound(k):
    K=k
    xsum = sum(arb(2)**i/arb(i).fac() for i in range(13, K//2+1))
    c013 = 1+14/SQ2
    Q13 = (a*c013/2).exp()*(-(14*PI*(K-13)).sqrt()).exp()
    rLm = (K*(K+1)/PI).sqrt()
    return 2*SQ2/14*(2*rLm*Q13 + (1/a)*(-(14*PI*(K-13)).sqrt()).exp())*xsum
for k in [99, 150, 400, 1000]:
    ex = mid2_exact(k); bd = mid2_bound(k)
    print(f"  k={k}: exact={ex.str(10)} bound(13)={bd.str(10)} exact/bound={(ex/bd).str(6)}")
    assert ex < bd

print("== (c) tail: exact box-sum vs (14) ==")
def tail_exact(k):
    return sum(arb(2)**i/arb(i).fac()*FH_box(k-i,i+1) for i in range(k//2+1, k+1))
def tail_bound(k):
    K=k
    xp = 2*(a/(2*SQ2)).exp()
    tail_sum1 = (a/(2*SQ2)).exp()*(arb(1).exp()*xp).exp()   # e^{a/(2 sqrt2)} e^{e x'}   <- (14)
    tail_sum2 = (2*arb(1).exp()).exp()                        # e^{2e}
    rLt = (K*(K+1)/PI).sqrt()
    return 2*SQ2/(K//2+1)*(2*rLt*tail_sum1 + tail_sum2/a)*(-arb(K)/2).exp()
def tail_bound_fixed(k):
    K=k
    xp = 2*(a/(2*SQ2)).exp()
    tail_sum1 = (a/2).exp()*(a/(2*SQ2)).exp()*(arb(1).exp()*xp).exp()  # WITH e^{a/2}
    tail_sum2 = (2*arb(1).exp()).exp()
    rLt = (K*(K+1)/PI).sqrt()
    return 2*SQ2/(K//2+1)*(2*rLt*tail_sum1 + tail_sum2/a)*(-arb(K)/2).exp()
for k in [99, 150, 400]:
    ex = tail_exact(k); bd = tail_bound(k); bf = tail_bound_fixed(k)
    print(f"  k={k}: exact={ex.str(10)}")
    print(f"        bound(14) as stated={bd.str(10)}  exact/bound={(ex/bd).str(6)}  VALID={ex<bd}")
    print(f"        bound(14)+e^(a/2)  ={bf.str(10)}  exact/bound={(ex/bf).str(6)}  VALID={ex<bf}")

print("== (d) per-term prefactor bounds ==")
for k in [99, 400, 10000]:
    # (k-i)(i+1) >= 14(k-13) for integer i in [13, k//2]
    bad = [i for i in range(13, k//2+1) if (k-i)*(i+1) < 14*(k-13)]
    rmax = max((2*arb(k-i)*(i+1)/PI).sqrt() for i in range(13, k//2+1))
    assert rmax <= (arb(k)*(k+1)/PI).sqrt()
    print(f"  k={k}: (k-i)(i+1)>=14(k-13) violations on [13,k//2]: {bad}; r_a<=sqrt(k(k+1)/pi) OK")
# corner: k=27 real-endpoint subtlety
print("  corner k=27: integer range is {13} only; (27-13)*14 =", 14*14, "= 14*(27-13): equality, OK")
print("  corner k=26..40 integer check:")
for k in range(26,41):
    bad = [i for i in range(13, k//2+1) if (k-i)*(i+1) < 14*(k-13)]
    if bad: print(f"    k={k}: VIOLATIONS at {bad}")
print("  (no output above = no integer violations)")
