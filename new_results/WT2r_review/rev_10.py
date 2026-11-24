
import sys
import sympy as sp
d, p = sp.symbols('d p', positive=True, integer=True)
I = sp.I; S2 = sp.sqrt(2)
tau1 = (-1-I)/2; tau2 = (-1+I)/2
def df(n): return sp.factorial(2*n-2)/(2**(n-1)*sp.factorial(n-1))
def Kfac(n):  # K_n without phase ambiguity, using i^{n+1} explicitly
    return (-1)**n * df(n) * I**(n+1)/(2*S2*sp.factorial(n))

# (A8) fixed: K_{2p} written with i^{2p+1} = i(-1)^p manually
K2p = df(2*p)*I*(-1)**p/(2*S2*sp.factorial(2*p))
nums = all(sp.simplify(Kfac(2*pp) - (df(2*pp)*I*(-1)**pp/(2*S2*sp.factorial(2*pp))))==0 for pp in range(1,40))
print("(A8) K_{2p} = i(-1)^p (4p-3)!!/(2 sqrt2 (2p)!) : verified p=1..39 exactly:", nums)
R = df(2*p)/(2*S2*sp.factorial(2*p))   # real positive
print("     K_{2p} = i(-1)^p R(p), R(p) real  =>  conj(rho_pp) = -rho_pp, rho_pp-conj = 2 rho_pp")

# (C5) fixed via gammasimp
ratio_lhs = sp.gammasimp(sp.binomial(sp.Rational(1,2), d+1)/sp.binomial(sp.Rational(1,2), d))
f = lambda n: (-1)**(n-1)*df(n)/(2**n*sp.factorial(n))
ratio_rhs = sp.simplify(f(d+1)/f(d))
print("(C5) binom(1/2,d) ratio:", ratio_lhs, "closed-form ratio:", ratio_rhs,
      "equal:", sp.simplify(ratio_lhs-ratio_rhs)==0, "; base d=1 equal:", sp.simplify(sp.binomial(sp.Rational(1,2),1)-f(1))==0)

print()
print("=== (I2) assembled: (2d-1) rho_pq = (p+1)tau1 rho_{p+1,q} + (q+1)tau2 rho_{p,q+1}, p>=q+2 ===")
q = d - p
# rho_pq = K_d (-1)^p 2 C(d,p) ; divide both sides by (-1)^p K_d * 2:
L = (2*d-1)*sp.binomial(d,p)
Kr = -I*(2*d-1)/(d+1)   # K_{d+1}/K_d from (A1)/(I1)
Rside = (p+1)*tau1*(-1)*Kr*sp.binomial(d+1,p+1) + (q+1)*tau2*Kr*sp.binomial(d+1,q+1)
# note rho_{p+1,q}/((-1)^p K_d 2) = -(K_{d+1}/K_d) C(d+1,p+1); rho_{p,q+1}/... = +(K_{d+1}/K_d) C(d+1,q+1)
Rside = sp.simplify(Rside)
diff = sp.simplify(L - Rside)
print("    [divided by (-1)^p K_d 2]  LHS-RHS =", diff)

print("=== (I3) assembled: (4p-3) rho_{p,p-1} = (p+1)tau1 rho_{p+1,p-1} + p tau2 (rho_pp - conj rho_pp) ===")
K2pm1 = Kfac(2*p-1)   # has i^{2p} = (-1)^p: rewrite
K2pm1 = df(2*p-1)*(-1)**(2*p-1)*I**(2*p)/(2*S2*sp.factorial(2*p-1))
# i^{2p} = (-1)^p :
K2pm1 = df(2*p-1)*(-1)*(-1)**p/(2*S2*sp.factorial(2*p-1))
rho_pm1 = K2pm1*(-1)**p*2*sp.binomial(2*p-1,p)          # K_{2p-1}(-1)^p 2 C(2p-1,p)
rho_p1m1 = K2p*(-1)**(p+1)*2*sp.binomial(2*p,p+1)       # K_{2p}(-1)^{p+1} 2 C(2p,p+1)
rho_pp = K2p*(-1)**p*sp.binomial(2*p,p)
L3 = (4*p-3)*rho_pm1
R3 = (p+1)*tau1*rho_p1m1 + p*tau2*(2*rho_pp)            # rho_pp - conj = 2 rho_pp (pure imaginary)
d3 = sp.simplify(sp.expand(L3 - R3))
# normalize the (-1)^p factors: divide by (-1)^p
d3n = sp.simplify(d3/(-1)**p)
print("    LHS-RHS (divided by (-1)^p) =", d3n)

print("=== (I4) assembled: Im hat_p = Im rho_pp, hat_p = (p+1)tau1 rho_{p+1,p}/(4p-1) ===")
K2p1 = df(2*p+1)*(-1)**(2*p+1)*I**(2*p+2)/(2*S2*sp.factorial(2*p+1))
# i^{2p+2} = (-1)^{p+1}
K2p1 = df(2*p+1)*(-1)*(-1)**(p+1)/(2*S2*sp.factorial(2*p+1))
rho_p1p = K2p1*(-1)**(p+1)*2*sp.binomial(2*p+1,p+1)
hat = (p+1)*tau1*rho_p1p/(4*p-1)
hat = sp.simplify(hat)
# extract imaginary part symbolically: write hat = X + iY with X,Y real (tau1 = (-1-i)/2)
hatX = sp.simplify((hat + hat.subs(I,-I))/2)   # wrong for non-real coeffs; do manual: tau1 -> (-1-I)/2 already explicit
# manual: expand hat into a + b*I
hate = sp.expand(hat)
hb = sp.simplify(hate.coeff(I,1)) if hate.coeff(I,1)!=0 else sp.simplify((hate - hate.subs(I,0))/I)
im_hat = sp.simplify((hate - hate.subs(I,0))/I)
im_rho = sp.simplify((sp.expand(rho_pp) - sp.expand(rho_pp).subs(I,0))/I)
print("    Im(hat_p) =", sp.simplify(im_hat))
print("    Im(rho_pp)=", sp.simplify(im_rho))
print("    equal:", sp.simplify(im_hat - im_rho)==0)
