# LAY1: checks of the reflection identity and the profile moments it yields (mpmath dps 60, sympy series).
import mpmath as mp, sympy as sp
mp.mp.dps = 60
g = lambda s: mp.log(1+mp.exp(-2*s))
Q = lambda f: mp.quad(f, [0, 1, 4, 16, 64, mp.inf])
z = mp.zeta
gam = {
 (0,1): Q(g), (1,1): Q(lambda s: s*g(s)), (0,2): Q(lambda s: g(s)**2)/2,
 (2,1): Q(lambda s: s**2*g(s))/2, (1,2): Q(lambda s: s*g(s)**2)/2, (0,3): Q(lambda s: g(s)**3)/6,
 (3,1): Q(lambda s: s**3*g(s))/6,
}
# (1) reflection identity G(l,m)+G(2m-l,m) = (1/2)[Gamma(-l/2)Gamma(l/2-m)/Gamma(-m) + 2/l + 2/(2m-l)] at a sample point
G = lambda l, m: Q(lambda s: mp.exp(l*s)*((1+mp.exp(-2*s))**m - 1))
l, m = mp.mpf('0.3'), mp.mpf('-0.45')
lhs = G(l, m) + G(2*m-l, m)
rhs = (mp.gamma(-l/2)*mp.gamma(l/2-m)/mp.gamma(-m) + 2/l + 2/(2*m-l))/2
print("reflection identity residual at (0.3,-0.45):", mp.nstr(lhs-rhs, 5))
# (2) series of the right side via sympy, coefficients of mu^k at lambda=0 and lambda^2 mu
L, M = sp.symbols('lambda mu')
R = (sp.gamma(-L/2)*sp.gamma(L/2-M)/sp.gamma(-M) + 2/L + 2/(2*M-L))/2
R0 = sp.simplify(sp.limit(R, L, 0))
print("R(0,mu) =", R0)
R0b = sp.simplify(R0.subs(sp.polygamma(0,-M), sp.polygamma(0,1-M)+1/M))
print("R(0,mu) rewritten:", R0b, " series:", sp.series(R0b, M, 0, 5))
# (3) moments vs closed forms
checks = [
 ("int_0^inf g        = zeta2/4 ", gam[(0,1)], z(2)/4),
 ("int_0^inf s g      = 3zeta3/16", gam[(1,1)], 3*z(3)/16),
 ("int_0^inf g^2/2    = zeta3/16", gam[(0,2)], z(3)/16),
 ("int_0^inf s^2 g/2  = 7zeta4/64", gam[(2,1)], 7*z(4)/64),
 ("int_0^inf s^3 g/6  = 7zeta5/32*? ", gam[(3,1)], 6*(1-mp.mpf(2)**-4)*z(5)/16/6),
 ("J = int(s g^2/2+g^3/6) = zeta4/32", gam[(1,2)]+gam[(0,3)], z(4)/32),
 ("deg2 reflection: g02+g11 = zeta3/4", gam[(0,2)]+gam[(1,1)], z(3)/4),
 ("deg3 reflection: 2g03+2g12+4g21 = zeta4/2", 2*gam[(0,3)]+2*gam[(1,2)]+4*gam[(2,1)], z(4)/2),
]
for name, a, b in checks:
    print(f"{name:45s} residual {mp.nstr(a-b, 5)}")
# (4) soft-max c4 and squircle c4 from the theorem's F_3
J = z(4)/32; ig3_soft = 2*(1/mp.mpf(3))*6*(1-mp.mpf(2)**-3)*z(4)/16; is2g = 2*gam[(2,1)]
print("softmax c4 = 4*(int_R g3 + 2J - int_R s^2 g/2) =", mp.nstr(4*(ig3_soft + 2*J - is2g), 30), " zeta4/4 =", mp.nstr(z(4)/4, 30))
ig3_sq = -2*(mp.mpf(4)/3)*6*(1-mp.mpf(2)**-3)*z(4)/16
c4sq = 4/mp.sqrt(2)*(ig3_sq + 2*J - is2g)
print("squircle c4 =", mp.nstr(c4sq, 30), " -33 sqrt2 zeta4/16 =", mp.nstr(-33*mp.sqrt(2)*z(4)/16, 30))
# (5) diagonal closed form G(t,t) = (1/4)[Gamma(-t/2)^2/Gamma(-t) + 4/t]
for t in [mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf('1.3')]:
    print("diagonal G(t,t) residual at t =", t, ":", mp.nstr(G(t, t) - (mp.gamma(-t/2)**2/mp.gamma(-t) + 4/t)/4, 5))
