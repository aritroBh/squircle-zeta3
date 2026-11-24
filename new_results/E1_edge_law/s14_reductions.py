# s14: exact symbolic reductions of the edge law Z_e = 6 al be (h^2/L)(3(al u + be v) - 2).
import sympy as sp
al, be, h, L, u, v = sp.symbols('alpha beta h L u v', positive=True)
Ze = 6*al*be*h**2/L*(3*(al*u + be*v) - 2)
# regular 2m-gon, K_1 = {sum_i |<n_i,x>| <= 1}: circumradius R = tan(x), x = pi/2m; h = R cos x, L = 2 R sin x,
# u = v = 1/2, al = be = w*|<n_0,V1>| = R sin(2x)
x, m = sp.symbols('x m', positive=True)
R = sp.tan(x)
reg = Ze.subs({al: R*sp.sin(2*x), be: R*sp.sin(2*x), h: R*sp.cos(x), L: 2*R*sp.sin(x), u: sp.Rational(1, 2), v: sp.Rational(1, 2)})
Cm = 48*m*sp.sin(x)**4*sp.cos(x)*(3*sp.sin(x)**2 - 1)
print("regular: 2m*Z_e - C_m =", sp.simplify(sp.expand_trig(2*m*reg - Cm)))
# rhombus |x/a| + |y/b| <= 1: al = be = 1, L = r, h = ab/r, u = a^2/r^2, v = b^2/r^2
a, b = sp.symbols('a b', positive=True); r = sp.sqrt(a**2 + b**2)
rh = Ze.subs({al: 1, be: 1, L: r, h: a*b/r, u: a**2/r**2, v: b**2/r**2})
print("rhombus: 4*Z_e - 24a^2b^2/r^3 =", sp.simplify(4*rh - 24*a**2*b**2/r**3))
# any parallelogram: al = be = 1 -> Z_e = 6 h^2/L (u + v = 1)
print("parallelogram edge: Z_e =", sp.simplify(Ze.subs({al: 1, be: 1, v: 1 - u})))
# affine-regular hexagon: al = be = 1/2 -> Z_e = -(3/4) h^2/L
print("affine-regular hexagon edge: Z_e =", sp.simplify(Ze.subs({al: sp.Rational(1, 2), be: sp.Rational(1, 2), v: 1 - u})))
# tangential polygon (inradius h): al = L sin(e0)/(2h), be = L sin(e1)/(2h), uL = h tan(e0/2), vL = h tan(e1/2)
e0, e1 = sp.symbols('epsilon0 epsilon1', positive=True)
Lt = h*(sp.tan(e0/2) + sp.tan(e1/2))
tg = Ze.subs({al: L*sp.sin(e0)/(2*h), be: L*sp.sin(e1)/(2*h), u: h*sp.tan(e0/2)/L, v: h*sp.tan(e1/2)/L})
F3 = sp.Rational(3, 2)*L*sp.sin(e0)*sp.sin(e1)*(3*sp.sin(e0/2)**2 + 3*sp.sin(e1/2)**2 - 2)
c0, c1 = sp.symbols('c0 c1', positive=True)   # c = tan(eps/2): sin eps = 2c/(1+c^2), sin^2(eps/2) = c^2/(1+c^2)
half = {sp.sin(e0): 2*c0/(1 + c0**2), sp.sin(e1): 2*c1/(1 + c1**2), sp.tan(e0/2): c0, sp.tan(e1/2): c1,
        sp.sin(e0/2)**2: c0**2/(1 + c0**2), sp.sin(e1/2)**2: c1**2/(1 + c1**2)}
print("tangential: Z_e - (3/2) L sin e0 sin e1 (3 sin^2(e0/2) + 3 sin^2(e1/2) - 2) =",
      sp.simplify(((tg - F3).subs(L, Lt)).subs(half)))
# large-m asymptotics of C_m
print("C_m series in 1/m:", sp.series(Cm.subs(x, sp.pi/(2*m)), m, sp.oo, 6))
import mpmath as mp
for mm in (10, 100, 1000):
    xx = mp.pi/(2*mm)
    print(f"m = {mm}: C_m * m^3 = {mp.nstr(48*mm*mp.sin(xx)**4*mp.cos(xx)*(3*mp.sin(xx)**2 - 1)*mm**3, 12)};  -3 pi^4 = {mp.nstr(-3*mp.pi**4, 12)}")
