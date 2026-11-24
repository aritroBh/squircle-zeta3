# P2N4R referee: (1) the a4 polynomial identity of Lemma g3:lem:a4; (2) an explicit C3 for Lemma g3:prop:n3 from Lemma g3:prop:n4.
import sympy as sp
from sympy import Rational as Q
out = []
def p(*s):
    print(*s); out.append(' '.join(map(str, s)))
Li2, S12, S22, w, l, z4 = sp.symbols('Li2 S12 S22 w ell zeta4')
u = w*l
f1 = -(1+u); f2 = 1+u-Li2; f3 = 2*w*S12 + (1+u)*(Li2-1)
f4 = Q(1,2)*Li2**2 - 2*S22 - u*Li2 - 2*w*S12 - Li2 + 1 + u
a4 = f4 - f1*f3 - Q(1,2)*f2**2 + f1**2*f2 - Q(1,4)*f1**4 + z4/2
p("a4(from f1..f4) - [zeta4/2 - 2 S22 + 2 w^2 ell S12 + (1-u^4)/4] =", sp.expand(a4 - (z4/2 - 2*S22 + 2*w**2*l*S12 + (1-u**4)/4)))
# check of the log-series coefficient formula itself
e = sp.symbols('e'); F1, F2, F3, F4 = sp.symbols('F1:5')
ls = sp.series(sp.log(1+F1*e+F2*e**2+F3*e**3+F4*e**4), e, 0, 5).removeO()
p("eps^4 coeff of log F - (f4 - f1 f3 - f2^2/2 + f1^2 f2 - f1^4/4) =", sp.expand(ls.coeff(e, 4) - (F4 - F1*F3 - F2**2/2 + F1**2*F2 - F1**4/4)))
# (2) C3: |E_n - (d-2 z a2 L)/n^3| <= |Q4|/n^4 + |R_n|, with |Q4|/z <= a2|(lam-3)^2-3| + 2|a3||lam-3| + 2|b4| + |p4|/z,
#     |a3| <= 387.5, |b4| <= 1937.5 (so |a2^2+2a4| = 2|b4|), 0<a2<2.15, |p4|/z = (1-z)|20z^2-2z+43|/6 <= 43/6 + ... ; lam in (L, L+1).
zz = sp.symbols('zz')
pz = (1-zz)*(20*zz**2-2*zz+43)/6
crit = [0, 1] + [r for r in sp.solve(sp.diff(pz, zz), zz) if r.is_real and 0 < r < 1]
pmax = max(pz.subs(zz, r) for r in crit)
p("max_{z in [0,1]} |p4|/z =", pmax, "=", float(pmax))
import mpmath as mm
L0 = 2.1972
def ratio(L):
    lam_hi = L + 1
    t1 = 2.15*max(3, (lam_hi-3)**2 - 3, 3 - min((x-3)**2 for x in (L, L+1, min(max(3, L), L+1))))
    return (t1 + 2*387.5*max(abs(L-3), abs(L+1-3)) + 2*1937.5 + float(pmax))/L**2
sup = max(ratio(L0 + k*0.001) for k in range(200000))
p("sup_{L>=log 9} |Q4|/(z L^2) <= %.1f (grid L in [2.1972, 202]; the ratio tends to 2.15 as L->oo)" % sup)
# |R_n| <= 1.2e5 z L^3/n^5 = 1.2e5 z (L/n) L^2/n^4, L/n <= log9/9 for n>=9
C3 = sup + 1.2e5*float(sp.log(9)/9)
p("explicit C3 <= %.0f  (so Lemma g3:prop:n3 holds with C3 = 3.1e4)" % C3)
open('ref_extra.log', 'w').write('\n'.join(out) + '\n')
