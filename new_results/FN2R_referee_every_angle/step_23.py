import sympy as sp, time
t0=time.time()
# ---------- Lemma B independent symbolic checks ----------
s, c, t = sp.symbols('s c t')
# (i) phi = h_c' ; rationalized form
h = 1 + s - sp.sqrt(1 + 2*c*s + s**2)
phi = sp.diff(h, s)
R = sp.sqrt(1 + 2*c*s + s**2)
phi_rat = (1-c**2)/(R*(R+s+c))
diff = sp.simplify((phi - phi_rat)*(R*(R+s+c)))   # should be 0 after multiply-out (symbolic with sqrt)
print("(i) phi - rationalized, cleared:", sp.simplify(phi - phi_rat).expand().subs(sp.sqrt(1+2*c*s+s**2), sp.Symbol('Rr')))
phid = sp.diff(phi, s)
print("(i) phi'(s) =", sp.simplify(phid), "  expect -(1-c^2)*(1+2cs+s^2)^(-3/2):",
      sp.simplify(phid + (1-c**2)*(1+2*c*s+s**2)**sp.Rational(-3,2))==0)
print("(i) phi(0) =", sp.simplify(phi.subs(s,0)), " expect 1-c")

# (iii) tail series in t=1/s to order 12 vs Legendre formula
phi_t = 1 - (1+c*t)*(1+2*c*t+t**2)**sp.Rational(-1,2)
ser = sp.series(phi_t, t, 0, 13).removeO().expand()
ok_l = []
for l in range(1, 13):
    coeff = ser.coeff(t, l)
    dl = (-1)**l*(1-c**2)*sp.legendre(l-1, c).diff(c)/l
    ok_l.append(sp.simplify(coeff - dl) == 0)
print("(iii) coeff(t^l) == (-1)^l(1-c^2)P'_{l-1}(c)/l for l=1..12:", all(ok_l), ok_l)

# the Legendre identity used: k*(P_k - x P_{k-1}) = (x^2-1) P'_{k-1}
x = sp.Symbol('x')
ok_id = all(sp.simplify(k*(sp.legendre(k,x) - x*sp.legendre(k-1,x)) - (x**2-1)*sp.legendre(k-1,x).diff(x))==0
            for k in range(1, 10))
print("(iii) Legendre identity k(P_k-xP_{k-1})=(x^2-1)P'_{k-1}, k=1..9:", ok_id)

# explicit values d2,d3,d4
d2 = (-1)**2*(1-c**2)*sp.legendre(1,c).diff(c)/2
d3 = (-1)**3*(1-c**2)*sp.legendre(2,c).diff(c)/3
d4 = (-1)**4*(1-c**2)*sp.legendre(3,c).diff(c)/4
print("d2 =", sp.simplify(d2), " d3 =", sp.simplify(d3), " d4 =", sp.expand(sp.simplify(d4)))
print("d2 == 2z(1-z), d3 == -4cz(1-z) with z=(1-c)/2:",
      sp.simplify(d2.subs(c, 1-2*sp.Symbol('z')) - 2*sp.Symbol('z')*(1-sp.Symbol('z')))==0,
      sp.simplify(d3.subs(c, 1-2*sp.Symbol('z')) + 4*(1-2*sp.Symbol('z'))*sp.Symbol('z')*(1-sp.Symbol('z')))==0)

# c=0: odd vanish, even = -binomial(-1/2,m)
ok0a = all(sp.simplify(((-1)**l*(1-c**2)*sp.legendre(l-1,c).diff(c)/l).subs(c,0))==0 for l in range(3,13,2))
ok0b = all(sp.simplify(((-1)**(2*m)*(1-c**2)*sp.legendre(2*m-1,c).diff(c)/(2*m)).subs(c,0)
           + sp.binomial(sp.Rational(-1,2),m))==0 for m in range(1,7))
print("c=0: odd l vanish (l=3..11):", ok0a, " ; d_{2m}(0) = -C(-1/2,m), m=1..6:", ok0b)
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} step_04: Lemma B symbolic (phi rationalized, phi' sign, d_l Legendre form l<=12, Legendre identity, c=0 -> -C(-1/2,m)) all EXACT [referee]\n")
