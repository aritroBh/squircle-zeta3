
import sys
sys.path.insert(0,'.')
import sympy as sp
v, c = sp.symbols('v c'); I = sp.I

# F(c/(v+i)^m) = c i^{-m} v^{m-2}/(v-i)^m  (exact rational algebra, m>=1)
# Then  F + c/(v-i)^m = c (i^{-m} v^{m-2} + 1)/(v-i)^m,  and numerator N(v) with N(i)=0
#  => pole order at v=i is <= m-1.  Hence top coeff is exactly -c.  Verify for m=1..80 exactly:
ok = True
for m in range(1, 81):
    R = c/(v+I)**m
    FR = sp.together(R.subs(v, 1/v)/v**2)
    closed = c * I**(-m) * v**(m-2)/(v-I)**m
    ok &= (sp.simplify(FR - closed) == 0)          # F formula exact
    N = I**(-m)*v**(m-2) + 1
    ok &= (sp.simplify(N.subs(v, I)) == 0)         # (v-i) divides N => order drops by >=1
print("F(c/(v+i)^m) = c i^{-m} v^{m-2}/(v-i)^m and N(i)=0 (top coeff exactly -c), m=1..80:", ok)

# ---- Lemma A ingredient 3 & 4 (rerun after crash) ----
sys.path.insert(0,'.')
import buildintegrand
from genintegrand import S_poly
okAB = all((lambda P: (not P[1]) and (not P[2]))(buildintegrand.AB_part_PF(k)) for k in range(1,7))
print("AB_part_PF(k) has empty +-i principal parts, k=1..6:", okAB)
okSQ = True
for d in range(2,10):
    for m in range(0,d+1):
        P = buildintegrand.sqrt_part_PF(d, S_poly(d,m))
        if len(P[1])>2*d-2 or len(P[2])>2*d-2: okSQ=False
print("sqrt_part_PF(d,m): +-i order <= 2d-2, d=2..9 all m:", okSQ)

# denominator structure of sqrt_part_PF for general d (from buildintegrand source):
# den = (v-a)(v+a)(v-b)(v+b) * v^{D-(2k-1)} (if D>2k-1) * (v^2+1)^{2k-2}; real poles simple.
# print for d=3, m=1 the PF to show structure:
from fastred import pf_to_sympy
P = buildintegrand.sqrt_part_PF(3, S_poly(3,1))
print("example sqrt_part (d=3,m=1): ipos length =", len(P[1]), " real poles:", [r.sympy() for r in P[0]])
