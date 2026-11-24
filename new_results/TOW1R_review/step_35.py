import sympy as sp
i=sp.I
p = sp.Symbol('p', integer=True, positive=True)
d = 2*p-1
def binom(a,b): return sp.binomial(a,b)
# rho_{pq}(0) proportional part: (-1)^P (1+[P>Q]) binom(D,P) * K_D. Ratio cancels K via K_{D+1}/K_D=-i(2D-1)/(D+1)
# We verify ratios using the closed form directly with a symbolic K ratio recurrence.
# rho(0;D,P) = K_D (-1)^P m_D(P) binom(D,P), m=1+[P>Q]. Build symbolic rho using factorials.
def K(D):  # K_D = (-1)^D (2D-3)!! i^{D+1}/(2 sqrt2 D!). Use ratio K_{D+1}/K_D = -i(2D-1)/(D+1); set K_d symbolic base.
    # represent K_D symbolically as Kd * prod; simpler: use explicit with double factorial
    from sympy import factorial, sqrt
    return (-1)**D * sp.factorial2(2*D-3) * i**(D+1) / (2*sqrt(2)*factorial(D))
def rho0(P,Q):
    D=P+Q; m = 2 if P>Q else 1
    return K(D)*(-1)**P*m*binom(D,P)
# ratio rA0 = rho_{p+1,p-1}/rho_{p,p-1}
rA0 = sp.simplify(rho0(p+1,p-1)/rho0(p,p-1))
print("rA0 =", sp.simplify(rA0), "  expected i*(2d-1)/(p+1) =", i*(2*d-1)/(p+1))
print("  match:", sp.simplify(rA0 - i*(2*d-1)/(p+1))==0)
rB0 = sp.simplify(rho0(p+1,p)/rho0(p,p-1))
print("rB0 =", sp.simplify(rB0), "  expected (2d-1)(2d+1)/((p+1)p) =", (2*d-1)*(2*d+1)/((p+1)*p))
print("  match:", sp.simplify(rB0 - (2*d-1)*(2*d+1)/((p+1)*p))==0)
# E3 ratio: rho_{p+1,p}/rho_{pp}  (diagonal rho_pp has m=1)
rB0_E3 = sp.simplify(rho0(p+1,p)/rho0(p,p))
print("E3 rB0 =", sp.simplify(rB0_E3), "  expected 2i(4p-1)/(p+1) =", 2*i*(4*p-1)/(p+1))
print("  match:", sp.simplify(rB0_E3 - 2*i*(4*p-1)/(p+1))==0)
# also verify conj(rho_pp(0)) = -rho_pp(0) (purely imaginary): conj K_{2p} relation
# K_{2p} = i^{2p+1}*real => conj = -K_{2p} iff i^{2p+1} pure imaginary => 2p+1 odd always true. So K_{2p} pure imaginary.
print("K_2p pure imaginary (p=3):", sp.re(K(6).expand(complex=True))==0, " K(6)=", K(6))
with open("<project>/research/08_theory/TOW1R_review/progress.log","a") as fh:
    fh.write("step_10: verified prove-script rho-ratios rA0,rB0 (E1,E2b) and E3 ratio + conj purity\n")
