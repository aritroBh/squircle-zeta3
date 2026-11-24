import sympy as sp, mpmath as mp
import formula
from formula import BLOCKS, B, r, q0, L2, V02
L = sp.sqrt(L2)
A0 = 1 - B*r; A1 = 1 - B; kap = r - 1
h2 = (L2*V02 - q0**2)/L2
L0 = -q0/L; L1 = L - L0
print('blocks parsed:', list(BLOCKS))
# 1. zeta(3) coefficient vs printed Z_e
Ze = 6*A0*A1*h2/L*(3*(A0*L0 + A1*L1)/L - 2)
print('c_z3/L^5 - Z_e =', sp.simplify(BLOCKS['zeta3']/L2**sp.Rational(5,2) - Ze))
# 2. Li3 coefficients vs printed Phi_A, Phi_L
gam = B*h2/L**3
PA = -6*gam*(3*A0*kap*q0 - L**2*(2*A0 + r*A1))
PL = -18*gam*(kap*(-A1**2 + A1*B*(2*kap-1) + 2*B**2*kap)*q0 + L**2*(A0*B + r*A1**2))
print('c_Li3(x)/L^5 - (2PhiA+PhiL) =', sp.simplify(BLOCKS['Li3(x)']/L2**sp.Rational(5,2) - (2*PA+PL)))
print('c_Li3(x/(x-1))/L^5 - PhiA =', sp.simplify(BLOCKS['Li3(x/(x-1))']/L2**sp.Rational(5,2) - PA))
# 3. Nielsen-Koelbig: 6A = 12Li3(x)+6Li3(x/(x-1))-6Li2(x)log(1-x)-log^3(1-x), A(x)=sum H_n x^n/n^2.
#    A'(x) = (Li2(x) + log(1-x)^2/2)/x.  Check derivative of RHS/6 symbolically, using Li2' = -log(1-x)/x,
#    Li3' = Li2/x and Landen Li2(x/(x-1)) = -Li2(x) - log(1-x)^2/2 (x<1).
x = sp.symbols('x'); Li2, Li3, Li3y = sp.symbols('Li2 Li3 Li3y'); lg = sp.log(1-x)
y = x/(x-1)
dLi2 = -lg/x; dLi3 = Li2/x
Li2y = -Li2 - lg**2/2               # Landen
dLi3y = Li2y/y*sp.diff(y, x)
dRHS = (12*dLi3 + 6*dLi3y - 6*(dLi2*lg + Li2*sp.diff(lg, x)) - 3*lg**2*sp.diff(lg, x))/6
print('d/dx[RHS/6] - (Li2+log^2(1-x)/2)/x =', sp.simplify(dRHS - (Li2 + lg**2/2)/x))
# both sides vanish at x=0 (Li_k(0)=0, log 1 = 0). Landen Li2 check numerically:
mp.mp.dps = 50
for xv in [mp.mpf('-3.5'), mp.mpf('-1'), mp.mpf('0.5'), mp.mpf('0.9')]:
    print('  Landen Li2 at', xv, mp.nstr(mp.polylog(2, xv/(xv-1)) + mp.polylog(2, xv) + mp.log(1-xv)**2/2, 3))
# numeric NK and special values
def Aser(xv):
    return mp.nsum(lambda n: mp.harmonic(n)*xv**n/n**2, [1, mp.inf])
def Arhs(xv):
    l = mp.log(1-xv)
    return (12*mp.polylog(3,xv) + 6*mp.polylog(3, xv/(xv-1)) - 6*mp.polylog(2,xv)*l - l**3)/6
for xv in [mp.mpf(-1), mp.mpf('0.5'), mp.mpf('-0.7'), mp.mpf('0.3')]:
    print('  NK at', xv, mp.nstr(Aser(xv) - Arhs(xv), 3))
z3 = mp.zeta(3); l2 = mp.log(2)
print('A(-1)+5/8 z3 =', mp.nstr(Aser(mp.mpf(-1)) + 5*z3/8, 3), ' A(1/2)-(z3-pi^2 log2/12) =', mp.nstr(Aser(mp.mpf(1)/2) - (z3 - mp.pi**2*l2/12), 3))
print('Li3(1/2) check', mp.nstr(mp.polylog(3, 0.5) - (7*z3/8 - mp.pi**2*l2/12 + l2**3/6), 3), ' Li3(-1)+3/4 z3', mp.nstr(mp.polylog(3,-1)+3*z3/4, 3))
# 4. the two other Euler-sum atoms of the engine
def Bser(xv): return mp.nsum(lambda n: mp.zeta(2) - mp.zeta(2, n+1) if False else sum(mp.mpf(1)/k**2 for k in range(1, int(n)+1))*xv**n/n, [1, mp.inf])
for xv in [mp.mpf('0.3'), mp.mpf('-0.6')]:
    l = -mp.log(1-xv)
    Bs = mp.nsum(lambda n: mp.harmonic(n, ) * 0 + mp.zeta(2) * 0 + sum(mp.mpf(1)/k**2 for k in range(1, int(n)+1))*xv**n/n, [1, mp.inf]) if False else None
mp.mp.dps = 40
for xv in [mp.mpf('0.3'), mp.mpf('-0.6'), mp.mpf('0.8')]:
    l = -mp.log(1-xv)
    s2 = 0; s1 = 0; Bs = 0; Cs = 0
    for n in range(1, 4000):
        s1 += mp.mpf(1)/n; s2 += mp.mpf(1)/n**2
        Bs += xv**n*s2/n; Cs += xv**n*s1**2/n
        if abs(xv)**n < mp.mpf(10)**(-45): break
    print('  sum H2_n x^n/n - (3Li3-2A+Li2 l) =', mp.nstr(Bs - (3*mp.polylog(3,xv) - 2*Aser(xv) + mp.polylog(2,xv)*l), 3),
          '   sum H_n^2 x^n/n - (Li3+Li2 l+l^3/3) =', mp.nstr(Cs - (mp.polylog(3,xv) + mp.polylog(2,xv)*l + l**3/3), 3))
# 5. orientation symmetry of the table formula: traverse an edge backwards
mp.mp.dps = 40
for (Bv, rv, P0, P1) in [(sp.Rational(3,8), 2, sp.Matrix([3,0]), sp.Matrix([sp.Rational(15,2), sp.Rational(3,4)])),
                          (sp.Rational(1,7), sp.Rational(5,2), sp.Matrix([1,2]), sp.Matrix([-1,3]))]:
    D = P1 - P0
    A0v = 1 - Bv*rv; A1v = 1 - Bv
    w0, w1, wB = 2, 3, 5
    ed = dict(B=Bv, r=rv, q0=P0.dot(D), L2=D.dot(D), V02=P0.dot(P0), A0=A0v, A1=A1v, wA0=w0, wA1=w1, wB=wB)
    Dr = P0 - P1
    edr = dict(B=Bv*rv, r=1/sp.nsimplify(rv), q0=P1.dot(Dr), L2=Dr.dot(Dr), V02=P1.dot(P1), A0=A1v, A1=A0v, wA0=w1, wA1=w0, wB=wB)
    print('  orientation: E(forward)-E(backward) =', mp.nstr(formula.E_formula(ed) - formula.E_formula(edr), 3))
