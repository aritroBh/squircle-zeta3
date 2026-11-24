
import sys
import sympy as sp
d, p = sp.symbols('d p', positive=True, integer=True)
I = sp.I; S2 = sp.sqrt(2)
tau1 = (-1-I)/2; tau2 = (-1+I)/2

def df(n):   # (2n-3)!! = (2n-2)!/(2^{n-1}(n-1)!)
    return sp.factorial(2*n-2)/(2**(n-1)*sp.factorial(n-1))
def K(n):    # K_n = (-1)^n (2n-3)!! i^{n+1}/(2 sqrt2 n!)
    return (-1)**n * df(n) * I**(n+1)/(2*S2*sp.factorial(n))

print("=== atomic identities (all symbolic) ===")
# (A1) K-ratio  (this is identity I1 itself)
r = sp.simplify(K(d+1)/K(d))
print("(A1)/(I1) K_{d+1}/K_d =", r, " == -i(2d-1)/(d+1):", sp.simplify(r - (-I*(2*d-1)/(d+1)))==0)
# (A2) absorption
print("(A2) (p+1)C(d+1,p+1) = (d+1)C(d,p):",
      sp.simplify((p+1)*sp.binomial(d+1,p+1) - (d+1)*sp.binomial(d,p))==0)
# (A3) symmetry
print("(A3) C(d,d-p) = C(d,p):",
      sp.simplify(sp.binomial(d,d-p).rewrite(sp.gamma) - sp.binomial(d,p).rewrite(sp.gamma))==0)
# (A4) tau facts
print("(A4) i(tau1-tau2)=1:", sp.simplify(I*(tau1-tau2))==1, " tau2-tau1=i:", sp.simplify(tau2-tau1)==I,
      " i*tau1=(1-i)/2:", sp.simplify(I*tau1-(1-I)/2)==0)
# (A5) duplication
print("(A5) C(2p,p) = 2C(2p-1,p):",
      sp.simplify(sp.binomial(2*p,p) - 2*sp.binomial(2*p-1,p))==0)
# (A6) absorption at (2p+1)
print("(A6) (p+1)C(2p+1,p+1) = (2p+1)C(2p,p):",
      sp.simplify((p+1)*sp.binomial(2*p+1,p+1) - (2*p+1)*sp.binomial(2*p,p))==0)
# (A7) C(2p,p+1) = C(2p,p) p/(p+1)
print("(A7) C(2p,p+1) = p/(p+1) C(2p,p):",
      sp.simplify(sp.binomial(2*p,p+1) - sp.Rational(1,1)*p/(p+1)*sp.binomial(2*p,p))==0)
# (A8) K_{2p} pure imaginary: K_{2p} = i*(-1)^p * R(p), R(p) = (4p-3)!!/(2 sqrt2 (2p)!)
R = df(2*p)/(2*S2*sp.factorial(2*p))
# i^{2p+1} = i*(-1)^p for integer p (standard: i^{2p} = (i^2)^p = (-1)^p)
check_i = all(sp.I**(2*pp+1) == sp.I*(-1)**pp for pp in range(1,60))
print("(A8) i^{2p+1} = i(-1)^p (p integer, checked p=1..59):", check_i,
      "; K_{2p} = i(-1)^p R(p):", sp.simplify(K(2*p) - I*(-1)**p*R)==0)
print("     => conj(K_{2p}) = -K_{2p} (pure imaginary), so rho_pp - conj rho_pp = 2 rho_pp")

# (C5) binom(1/2,d) = (-1)^{d-1} (2d-3)!!/(2^d d!) via ratio recurrence + base case
ratio_lhs = sp.simplify(sp.binomial(sp.Rational(1,2), d+1)/sp.binomial(sp.Rational(1,2), d))
f = lambda n: (-1)**(n-1)*df(n)/(2**n*sp.factorial(n))
ratio_rhs = sp.simplify(f(d+1)/f(d))
print("(C5) ratio binom(1/2,d+1)/binom(1/2,d) =", ratio_lhs, "; closed-form ratio =", ratio_rhs,
      "; equal:", sp.simplify(ratio_lhs-ratio_rhs)==0,
      "; base d=1:", sp.simplify(sp.binomial(sp.Rational(1,2),1) - f(1))==0)
