import mpmath
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, nstr, catalan
mp.dps = 70
# fit (P-8)*p = c1 + c2/p + c3/p^2 from p=12,16,20 (real parts)
vals = [(mpf(12), mpf('1.11875104547822741838032431427')),
        (mpf(16), mpf('1.11875104671492038555794907271')),
        (mpf(20), mpf('1.11875104671543652916815781961'))]
from mpmath import matrix, lu_solve
A = matrix(3,3); b = matrix(3,1)
for i,(p,v) in enumerate(vals):
    A[i,0], A[i,1], A[i,2] = 1, 1/p, 1/p**2
    b[i] = v
c = lu_solve(A,b)
print("c1 =", nstr(c[0],45))
print("c2 =", nstr(c[1],10), " c3 =", nstr(c[2],10))
print("c1_pred = 1.118751046715436736371392659321737046796877152219")
# I_c = int cos u ln cos u over [0, pi/4]
Ic = mpmath.quad(lambda u: cos(u)*log(cos(u)), [0, pi/8, pi/4])
print("I_c =", nstr(Ic,55))
# PSLQ Ic
L2 = log(2)
lasinh1 = log(1+sqrt(2))
cand = {'1':mpf(1),'C':catalan,'pi':pi,'L2':L2,'piL2':pi*L2,'asinh1':lasinh1,'L2^2':L2**2,'pi^2':pi**2}
basis = list(cand.values()); names = list(cand.keys())
rel = mpmath.pslq(basis+[Ic], maxcoeff=10**8)
print("PSLQ I_c:", rel)
if rel and rel[-1]!=0:
    print(" + ".join(f"({r})*{nm}" for r,nm in zip(rel,names) if r) + f" + ({rel[-1]})*Ic")
# c1 = 12 ln2 - 8 - 8 sqrt2 Ic check
c1 = mpf('1.118751046715436736371392659321737046796877152219')
print("12ln2-8-8sqrt2 Ic =", nstr(12*L2-8-8*sqrt(2)*Ic,45))
