"""Independent derivation of varpi2 (and control varpi1) by the classical second-order Laplace formula
in derivative form (not Gaussian moments), in the real variable r:
  int g e^{-Lam f} dr ~ e^{-Lam f0} sqrt(2pi/(Lam f2)) [ g + (1/Lam)( g''/(2f2) - g' f3/(2 f2^2)
                                                         - g f4/(8 f2^2) + 5 g f3^2/(24 f2^3) ) ],
with f = -psi, f2=f'', f3=f''', f4=f'''' at r=1 (principal root; Re f2 > 0).
q=2: s = 2 lam r w, lam = sqrt(N/pi), Lam = pi lam; integrand lam^{-3/2} p(r)(1+eps) e^{Lam psi},
     p = r^{-5/2} e^{-i pi/(2 r^2)} (constant dropped), eps = (1/lam)[ -3/(8 r w) + 1/(2 r w) - pi/(3 r^3 w^3) ] + O(lam^-2)
     (Stirling Gamma(z)/Gamma(z+3/2) = z^{-3/2}(1 - 3/(8z) + ...), z = s/2;  1 + x/2, x = 2/s;  e^{-N x^3/3}).
     relative correction = C/Lam  ->  varpi2 = C/sqrt(pi)   (1/Lam = 1/(sqrt(pi) sqrt(N))).
q=1: s = lam1 r w, lam1 = sqrt(2N/pi), Lam1 = sqrt(pi N/2) = (pi/2) lam1; p = r^{-5/2} e^{-i pi/(4 r^2)},
     eps = (1/lam1)[ -3/(4 r w) - pi/(6 r^3 w^3) ];  varpi1 = C1 * sqrt(2/pi)."""
import sympy as sp
r = sp.symbols('r', positive=True)
w = sp.exp(sp.I*sp.pi/4)
psi = sp.I*w*r - sp.conjugate(w)/r
f = -psi
def laplace_rel(g, f):
    d = lambda e, n: sp.diff(e, r, n).subs(r, 1)
    f2, f3, f4 = d(f, 2), d(f, 3), d(f, 4)
    g0, g1, g2 = d(g, 0), d(g, 1), d(g, 2)
    assert sp.simplify(d(f, 1)) == 0
    return sp.simplify((g2/(2*f2) - g1*f3/(2*f2**2) - g0*f4/(8*f2**2) + 5*g0*f3**2/(24*f2**3))/g0), f2
# q = 2
p2 = r**sp.Rational(-5, 2)*sp.exp(-sp.I*sp.pi/(2*r**2))
relL, f2 = laplace_rel(p2, f)
eps2 = (-sp.Rational(3, 8)/(r*w) + sp.Rational(1, 2)/(r*w) - sp.pi/(3*r**3*w**3)).subs(r, 1)
C2 = relL + sp.pi*eps2                       # 1/lam = pi/Lam
v2 = sp.nsimplify(sp.expand_complex(C2/sp.sqrt(sp.pi)))
v2n = sp.N(C2/sp.sqrt(sp.pi), 40)
rsts = sp.sqrt(2)*((12+33*sp.pi+2*sp.pi**2) + sp.I*(12-33*sp.pi+2*sp.pi**2))/(48*sp.sqrt(sp.pi))
print("f''(1) =", sp.simplify(f2), "  (RSTS: -psi''(1) = 2 conj(w))")
print("varpi2 (referee) =", v2n)
print("varpi2 (RSTS)    =", sp.N(rsts, 40))
print("difference       =", sp.N(C2/sp.sqrt(sp.pi) - rsts, 50))
# leading constant check: lam^{-3/2} * (-(sqrt(pi)/2) w^{-3/2}) * p(1) e^{Lam psi(1)} sqrt(2 pi/(Lam f2))
lead = -(sp.sqrt(sp.pi)/2)*sp.exp(-3*sp.I*sp.pi/8)*sp.exp(-sp.I*sp.pi/2)*sp.sqrt(2*sp.pi/f2)
print("leading phase/amplitude factor (x Lam^{-1/2} lam^{-3/2}):", sp.N(lead, 20), " expected -(pi/2) e^{-3 i pi/4} =", sp.N(-(sp.pi/2)*sp.exp(-3*sp.I*sp.pi/4), 20))
print("psi(1) =", sp.simplify(psi.subs(r, 1)))
# q = 1 control
p1 = r**sp.Rational(-5, 2)*sp.exp(-sp.I*sp.pi/(4*r**2))
relL1, _ = laplace_rel(p1, f)
eps1 = (-sp.Rational(3, 4)/(r*w) - sp.pi/(6*r**3*w**3)).subs(r, 1)
C1 = relL1 + sp.pi/2*eps1                    # 1/lam1 = (pi/2)/Lam1
v1 = C1*sp.sqrt(2)/sp.sqrt(sp.pi)            # 1/Lam1 = sqrt(2)/(sqrt(pi) sqrt(N))
paper = ((24+12*sp.pi+sp.pi**2) + sp.I*(24-12*sp.pi+sp.pi**2))/(48*sp.sqrt(sp.pi))
print("control varpi1 (referee) =", sp.N(v1, 30), "  paper:", sp.N(paper, 30), "  diff:", sp.N(v1-paper, 40))
