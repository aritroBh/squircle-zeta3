"""RST1 step_18: certified validation of Lemmas R2 (Gaussian F_H) and R3 (box, rate lam).
Exact F_H^{lam}(j,q) by arb quadrature vs the analytic bounds; prints ratios.
Also validates the convexity lemma (Q_i <= Q_2 on [2,k/2]; Q_13 >= Q_{k/2}).
"""
from flint import arb, acb, ctx
ctx.dps = 50
PI = arb.pi(); SQ2 = arb(2).sqrt(); sqpi = PI.sqrt()
a = PI/(2*SQ2)
W = acb(1,1)/SQ2

def FH_exact(j, q, lam):
    """arb quadrature of int_0^inf e^{-lam r} |H(s)| |s|^j/|s+q|^{j+1} dr, s=1/2+r W."""
    jj = arb(j); qq = arb(q)
    def f(r, analytic):
        s = arb(1)/2 + r*W
        H = sqpi/4*(s/2).gamma()*(s/2+arb(3)/2).rgamma()
        return (-lam*r).exp()*abs(H)*abs(s)**jj/abs(s+qq)**(jj+1)
    value = acb(0)
    pts=[0,1,4,16,64,200]
    for lo,hi in zip(pts,pts[1:]):
        piece = acb.integral(f, lo, hi, abs_tol=arb("1e-70"), rel_tol=arb("1e-60"), eval_limit=80000)
        assert piece.is_finite()
        value += piece
    # tail allowance: |H|<=sqrt(pi/2) x^-3/2, (|s|/|s+q|)^j<=1, 1/|s+q|<=sqrt2/q
    xR = arb(1)/2 + 200/SQ2
    tail = sqpi/arb(2).sqrt()**-1 * 0 + sqpi/SQ2 * xR**arb('-1.5')*SQ2/qq*(-lam*200).exp()/lam
    return value.real + arb(0, tail.upper())

def FH_gauss(j, q, lam):
    j,q = arb(j),arb(q)
    rL=(j*q/(SQ2*lam)).sqrt(); c0=1+q/SQ2; Sp=(2*SQ2*lam*j*q).sqrt()
    assert rL>c0
    CL=2**arb('1.5')/(3*lam); c2=8*SQ2*j*q/(4*rL+c0)**3
    mid=((rL-c0)/(2*SQ2))**arb('-1.5')*(2*PI/c2).sqrt()
    return sqpi/q*(CL*(-Sp/4).exp()+mid+(SQ2*rL)**arb('-1.5')/lam)*(lam*c0/2).exp()*(-Sp).exp()

def FH_box(j,q,lam):
    j,q=arb(j),arb(q)
    rL=(j*q/(SQ2*lam)).sqrt(); c0=1+q/SQ2; Sp=(2*SQ2*lam*j*q).sqrt()
    return 2*SQ2/q*(2*rL*(lam*c0/2).exp()+1/lam)*(-Sp).exp()

print("== Lemma R2 (Gaussian) and R3 (box) validation: exact/bound ==")
cases = [(99,2,a),(98,2,a),(400,2,a),(87,13,a),(97,3,a),(91,9,a),(99,1,3*a),(400,1,3*a),(50,5,a),(10,2,a),(200,6,a)]
worst=0
for (j,q,lam) in cases:
    Fe=FH_exact(j,q,lam); Fg=FH_gauss(j,q,lam); Fb=FH_box(j,q,lam)
    rg=(Fe/Fg); rb=(Fe/Fb)
    worst=max(worst, rg.upper(), rb.upper())
    print(f" j={j} q={q} lam/a={mp if False else (lam/a).str(6)}: exact={Fe.str(8)} gauss/bound ratio={rg.str(6)} box ratio={rb.str(6)}", flush=True)
print(f"worst exact/bound ratio: {worst.str(8)} (must be <1)", flush=True)

print("== convexity lemma ==")
# Q_i = e^{a(i+1)/(2 sqrt2)} e^{-sqrt(pi (k-i)(i+1))}; check Q_i<=Q_2 for i in [2,k/2], k=99..400
def Q(i,k):
    i,k=arb(i),arb(k)
    return (a*(i+1)/(2*SQ2)).exp()*(-(PI*(k-i)*(i+1)).sqrt()).exp()
bad=0
for k in [99,150,300,400]:
    Q2=Q(2,k)
    for i in range(2,k//2+1):
        if Q(i,k)>Q2: bad+=1
print(f"violations Q_i<=Q_2 (i in [2,k/2]): {bad}", flush=True)
# Q_13 >= Q_{k/2}
bad2=0
for k in range(26,401):
    if Q(k//2,k) > Q(13,k): bad2+=1
print(f"violations Q_13>=Q_(k/2) for k in [26,400]: {bad2}", flush=True)
