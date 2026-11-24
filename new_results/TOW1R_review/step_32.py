import sympy as sp
v = sp.symbols('v'); i=sp.I; sq2=sp.sqrt(2)
a=1+sq2; b=sq2-1
DL1 = 1/(v-b)+1/(v+a)-1/v
DL2 = 1/(v-a)+1/(v+b)-1/v
# key fold lemma: DL1(1/v)/v^2 = -DL2(v)
lhs = sp.simplify(DL1.subs(v,1/v)/v**2)
rhs = sp.simplify(-DL2)
print("DL1(1/v)/v^2 =", sp.simplify(lhs))
print("-DL2(v)      =", sp.simplify(rhs))
print("fold lemma DL1(1/v)/v^2 == -DL2(v):", sp.simplify(lhs-rhs)==0)

# E1 non-adjacent symbolic (d,p independent symbols, q=d-p, gap>=2)
d,p = sp.symbols('d p')
q = d-p
tau10=(-1-i)/2; tau20=(-1+i)/2; tau11=sp.Rational(3,4)*(1-i); tau21=-sp.Rational(3,4)*(1+i)
rA0 = i*(2*d-1)/(p+1)      # rho_{p+1,q}(0)/rho_{p,q}(0)
rB0 = -i*(2*d-1)/(q+1)     # rho_{p,q+1}(0)/rho_{p,q}(0)
LHS = -i*(d-2)
RHS = ( (p+1)*( tau11*rA0/(2*d-1) + tau10*(-i*(d-1))*rA0/(2*d-2) )
      + (q+1)*( tau21*rB0/(2*d-1) + tau20*(-i*(d-1))*rB0/(2*d-2) ) )
print("E1 non-adjacent remainder:", sp.simplify(sp.expand(LHS-RHS)))
with open("<project>/research/08_theory/TOW1R_review/progress.log","a") as fh:
    fh.write("step_07: fold lemma DL1(1/v)/v^2=-DL2(v) PROVED; E1 non-adjacent re-derived symbolically -> 0\n")
