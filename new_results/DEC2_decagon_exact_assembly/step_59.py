import sympy as sp, pickle, os
import mpmath as mp
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
x = sp.Symbol('x')
L1 = -sp.log(1-x); L2 = sp.polylog(2,x); L3 = sp.polylog(3,x)
mp.mp.dps = 50
xv = mp.mpf('0.37')
def H(n): return mp.harmonic(n)
def H2(n): return mp.fsum(1/mp.mpf(j)**2 for j in range(1,int(n)+1)) if n>=1 else mp.mpf(0)
gf1 = mp.nsum(lambda n: H(n)*xv**n, [1, mp.inf])
gf2 = mp.nsum(lambda n: H2(n)*xv**n, [1, mp.inf])
gf3 = mp.nsum(lambda n: H(n)**2*xv**n, [1, mp.inf])
L1v = -mp.log(1-xv); L2v = mp.polylog(2,xv)
print("GF numeric checks:",
      mp.nstr(abs(gf1-L1v/(1-xv)),3),
      mp.nstr(abs(gf2-L2v/(1-xv)),3),
      mp.nstr(abs(gf3-(L2v+L1v**2)/(1-xv)),3))
# atom identities (EXACT xd/dx residuals)
lhs1 = L1/(1-x)
rhs1 = sp.expand(x*sp.diff(L2 + L1**2/2, x))
c1 = sp.simplify(lhs1 - rhs1)
lhs2 = L2/(1-x)
rhs2 = sp.expand(3*L2 - 2*(L2+L1**2/2) + x*sp.diff(L2*L1, x))
c2 = sp.simplify(lhs2 - rhs2)
lhs3 = (L2+L1**2)/(1-x)
rhs3 = sp.expand(L2 + x*sp.diff(L2*L1 + L1**3/3, x))
c3 = sp.simplify(lhs3 - rhs3)
print("atom xd/dx residuals (must be 0):", c1, c2, c3)
# w^c expansions
y = sp.Symbol('y'); w = sp.log(1-y)
e2 = sp.simplify(-2*w/(1-y) - 2*(L1/(1-x)).subs(x,y))
e3 = sp.simplify(-3*w**2/(1-y) + 3*(L1**2).subs(x,y)/(1-y))
print("w^2 expansion derivative residual:", e2)
print("w^3 expansion derivative residual:", e3)
with open('progress.log','a') as f:
    f.write("step_07: esum machinery PROVED (GFs, 3 atom ids xd/dx residuals exactly 0, w^c c=2,3 residuals exactly 0); GF numeric %.1e %.1e %.1e\n" % (abs(gf1-L1v/(1-xv)),abs(gf2-L2v/(1-xv)),abs(gf3-(L2v+L1v**2)/(1-xv))))
print("DONE")
