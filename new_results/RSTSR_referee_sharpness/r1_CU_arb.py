"""RSTSR referee: independent Arb recomputation of C_U (RSTS (ii)), written from the referee's own
term-by-term derivation (see REFEREE.md sec. 3).  Bounds used (CX6, refereed):
 (13) on Gamma: |H|<=2, |K|<=2e^{-a r}, |K-E|<=2e^{-3a r}, a=pi/(2 sqrt2)
 (14) F(j,q;lam) <= (4/lam) e^{lam/2 + lam q/sqrt2} e^{-sqrt(3 lam j q/sqrt2)}
 (16) sum_{i>k/2} x^i/i! <= e^{e x} e^{-k/2}
 (17) (2/pi)|I3| <= T0 e^{-sqrt(9 pi k/4)}
Each term t_j(N) <= c_j * e^{-(3/2)sqrt(pi(N-3))} for N>=100; C_U = sum c_j.
Also checks the exponent comparisons on a grid AND by the referee's hand argument (printed)."""
from flint import arb, ctx
ctx.prec = 256
pi = arb.pi(); s2 = arb(2).sqrt(); a = pi/(2*s2); e = arb(1).exp()
def pref(q, lam):           # prefactor in (14)
    return 4/lam*(lam/2 + lam*q/s2).exp()
x = 2*(pi/4).exp()
twopi = 2/pi
# I1 = int_Gamma H (K-E) A2 : |H||K-E||A2| <= 2*2e^{-3ar}*2 rho^{k-1}/|s+2|^k  -> 8 F(k-1,2;3a)
cI1 = twopi*8*pref(2, 3*a)
# I2 = int_Gamma H K A3 : |H K| <= 4 e^{-ar};  (B3): 2F(k,3)+4F(k-1,3)+2 sum_{i>=2} 2^i/i! F(k-i,i+1)
#   pref(i+1,a) = (4/a) e^{a/2} e^{pi (i+1)/4};  sum_{i>=2} (2^i/i!) e^{pi i/4} <= e^{x} - 1 - x  (use e^x)
pa = 4/a*(a/2).exp()*(pi/4).exp()
cI2_lead = twopi*4*(2*pref(3, a) + 4*pref(3, a))
cI2_mid  = twopi*4*2*pa*x.exp()
cI2_tail = twopi*4*2*pa*(e*x).exp()
T0 = 32/(3*pi*a)*(3*a/2 + 3*pi/4).exp()
# elementary: |J_k| <= e^{-(N-1)/2};  (2/pi)(seg_Q + seg_A2) <= e^{-(N-1)/2}
cel = arb(2)
CU = cI1 + cI2_lead + cI2_mid + cI2_tail + T0 + cel
for nm, v in [("I1", cI1), ("I2 lead (i=0,1)", cI2_lead), ("I2 mid (2<=i<=k/2)", cI2_mid),
              ("I2 tail (i>k/2)", cI2_tail), ("I3 = T0", T0), ("elementary", cel), ("C_U", CU)]:
    print(f"{nm:22s} {v.str(18, radius=True)}")
print("C_U upper bound:", CU.upper().str(15), "  < 1.0584e7 ?", bool(CU < arb("1.0584e7")), "  < 1.06e7 ?", bool(CU < arb("1.06e7")))
# sharper variant: sum_{i>=2} 2^i/i! e^{pi i/4} = e^x - 1 - x exactly
cI2_mid_sharp = twopi*4*2*pa*(x.exp() - 1 - x)
print("C_U with e^x-1-x       ", (CU - cI2_mid + cI2_mid_sharp).str(15))
# exponent checks used: each exponent >= E(N) = 1.5 sqrt(pi (N-3))
import mpmath as mp
mp.mp.dps = 30
bad = []
for N in list(range(100, 5001)) + [10**j for j in range(4, 13)]:
    k = N - 1; E = 1.5*mp.sqrt(mp.pi*(N-3))
    exps = {"I1": mp.sqrt(9*mp.pi*(N-2)/2), "F(k,3)": mp.sqrt(9*mp.pi*k/4), "F(k-1,3)": mp.sqrt(9*mp.pi*(k-1)/4),
            "I3": mp.sqrt(9*mp.pi*k/4), "k/2": mp.mpf(k)/2, "(N-1)/2": mp.mpf(N-1)/2}
    # mid-range: min over 2<=i<=k/2 of 3pi(k-i)(i+1)/4
    mn = min((k-i)*(i+1) for i in (2, k//2)) if k >= 4 else None
    exps["mid min"] = mp.sqrt(3*mp.pi*mn/4)
    for nm, v in exps.items():
        if v < E - mp.mpf(10)**-20: bad.append((N, nm))
print("exponent comparison failures on N in [100,5000] U {1e4..1e12}:", bad[:5], len(bad))
# identity (k-i)(i+1)-3(k-2) = (i-2)(k-3-i)
import sympy as sp
i_, k_ = sp.symbols('i k')
print("identity:", sp.expand((k_-i_)*(i_+1) - 3*(k_-2) - (i_-2)*(k_-3-i_)) == 0)
# elementary terms numerically at N=100: J_k and the two segments vs e^{-(N-1)/2}
for N in (100, 101, 200):
    k = N-1
    J = mp.log(2)**N/mp.factorial(N)
    seg = 2/(mp.pi*k)*mp.mpf(3)**(-k) + 6/mp.pi*mp.mpf(2)**(-2*k)/(k-1)
    print(f"N={N}: J={mp.nstr(J,5)} seg={mp.nstr(seg,5)} e^-(N-1)/2={mp.nstr(mp.exp(-(N-1)/mp.mpf(2)),5)} E-bound={mp.nstr(mp.exp(-1.5*mp.sqrt(mp.pi*(N-3))),5)}")
