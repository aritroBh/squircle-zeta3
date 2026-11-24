with open('progress.log','a') as f:
    f.write("turn12-14 R4 (H2): x^2|g(x)| max on (0,1] = 0.999998915112 (200000 pts), ->1 as x->0 (tiny-x max 1.0): <=1 CONFIRMED.\n"
            "  c_k>0 (k<=60); sum c_k = 1 EXACT (sympy); partial sums track 1-1/sqrt(pi K); nsum's 0.9985 was an nsum artifact.\n"
            "  series identity resid ->0 with truncation (1e-9 @K=2e5); beta-integral identity EXACT by B-function def (resid 2.3e-22 quad).\n"
            "  min 2Re(w_k)/|w_k| = 1.4142135 >= 1 (dense x in (0,1], k up to 1e8). I0/bound ratio ~0.27 at N=100,400,1e3,1e4.\n"
            "  e^RePhi = e^{-ax} x^N (1+rt2 x+x^2)^{-N/2} resid 4.6e-48. Proof chain re-derived by hand: no gap (Weierstrass M, Jensen, |w|<=1, (x-1)^2).\n"
            "  R4 VERDICT: HOLDS (proof complete + dense numerics).\n")
# R5(a): reproduce certify_H3 + subdivision + dense truth
from flint import arb, ctx
ctx.dps = 60
rt2 = arb(2).sqrt()
def iv(lo,hi):
    lo=arb(lo).lower(); hi=arb(hi).upper(); return arb((lo+hi)/2,(hi-lo)/2)
def Hp3(Tb):
    den = 1+rt2*Tb+Tb**2
    q = ((1+rt2*Tb)/den).sqrt()
    w = 1/(den*(1+q))
    Gp = (rt2/3)*(rt2+2*Tb)/den**2
    hp = rt2*w**2*Tb + rt2*w*Tb*(2+rt2*Tb)/(q*den**2)
    return Gp - hp
v = Hp3(iv(0.0,0.13))
print(f"R5(a) single box Hp3([0,0.13]) = {v}  lower = {float(arb(v.lower()))}")
mn = arb(9)
for i in range(1300):
    a0 = arb(i)/10000; b0 = arb(i+1)/10000
    vv = Hp3(iv(a0,b0))
    mn = min(mn, arb(vv.lower()))
print(f"      1300 sub-boxes of [0,0.13]: min of lower bounds = {float(mn)}")
