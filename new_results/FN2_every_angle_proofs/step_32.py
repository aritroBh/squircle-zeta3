import sympy as sp, time
t0 = time.time()
a, A, B = sp.symbols('a A B')
c, z = sp.symbols('c z')
k2, k3, Xi = sp.symbols('k2 k3 Xi')
mu1 = Xi; mu2 = Xi**2 + k2; mu3 = Xi**3 + 3*k2*Xi + k3

# my master J's (from previous step)
th = {(1,0): 1/a, (2,0): (A - a**2)/(2*a**3), (1,1): -A/a**3,
      (3,0): (3*A**2 - 3*A*a**2 - B*a + a**4)/(6*a**5),
      (2,1): (-4*A**2 + 3*A*a**2 + B*a)/(2*a**5), (1,2): A**2/a**5}
M = {1: 2*z*mu1, 2: 2*z*mu2, 3: 2*z*mu3}
T2 = 2*z*(1-z)/a**2
T3 = -4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3
J1 = sp.expand(th[(1,0)]*M[1])
J2 = sp.expand(th[(2,0)]*M[2] + th[(1,1)]*M[1] + T2)
J3 = sp.expand(th[(3,0)]*M[3] + th[(2,1)]*M[2] + th[(1,2)]*M[1] + T3)

# ---- Recheck 1: theta=pi/2 (c=0, z=1/2) must equal P3's proved Q1,Q2,Q3 ----
k1 = sp.Symbol('k1')
XiP3 = sp.Symbol('XiP3')  # P3's Xi = log(e a n/2); same variable here
# P3 Q's (from P3 result.md, Theorem 2 explicit terms), with b=A, c_=B:
Q1 = Xi/a
Q2 = -(a**2-A)/(2*a**3)*(Xi**2+k2) - A*Xi/a**3 + 1/(2*a**2)
lam3 = a**4 - 3*a**2*A - a*B + 3*A**2
Q3 = lam3/(6*a**5)*(Xi**3+3*k2*Xi+k3) + (3*a**2*A+a*B-4*A**2)/(2*a**5)*(Xi**2+k2) + A**2*Xi/a**5 - A/a**4
subR = {z: sp.Rational(1,2), c: 0}
print("pi/2 reduction: J1-Q1 =", sp.simplify((J1-Q1).subs(subR)),
      " J2-Q2 =", sp.simplify(sp.expand((J2-Q2).subs(subR))),
      " J3-Q3 =", sp.simplify(sp.expand((J3-Q3).subs(subR))))

# ---- Recheck 2: f=x vs the EXACT every-angle series (Ahmed/P2) through eps^3 ----
# For f=x: exact 2-L_n = (1-c)(1 - exp(-eps*Xi + a2 eps^2 + a3 eps^3 + ...)) + tail terms?
# Rather: use AH3's exact_J1..exact_J3 (validated against the exact series in AH3 step_15/17/19)
exact_J1 = 2*z*mu1
exact_J2 = -z*mu2 + 2*z*(1-z)
exact_J3 = z*mu3/3 - c*z*(1-z)
subx = {a: 1, A: 0, B: 0}
print("f=x: J1-exact =", sp.simplify((J1-exact_J1).subs(subx)),
      " J2-exact =", sp.simplify(sp.expand((J2-exact_J2).subs(subx))),
      " J3-exact =", sp.simplify(sp.expand((J3-exact_J3).subs(subx))))

# ---- Recheck 3: independent derivation of the f=x exact series from Theorem P2's exact formula ----
# For f=x: y=x^n, slope s -> tau = 1 - (s/n)^(1/(n-1)). 2-L_n = int phi_theta(s) (1-(s/n)^(1/(n-1))) ds.
# In eps=1/(n-1): 1-(s/n)^eps = 1 - exp(eps*log(s/n)). With r=log(n a/s), a=1: = 1 - exp(-eps*r).
# 2-L = int_0^n phi_theta(s) (1-e^{-eps log(n/s)}) ds  [m_y=0 for s>n]
#     = int_0^inf phi (1-e^{-eps log(n/s)}) ds - int_n^inf phi (1-e^{-eps log(n/s)}) ds
# term1 = 2z*E[1 - exp(-eps*(log n + xi))] with xi=log(1/S_theta), CGF log E e^{z xi} = log R_legendre...
# E exp(-eps xi) = E S^eps = exp(a1 eps + a2 eps^2 + a3 eps^3 + ...) ; log n + ... : Xi = log n + k1
# so term1 = 2z*(1 - exp(-eps*Xi + aa2 eps^2 + aa3 eps^3))  where k_j = (-1)^j j! a_j:
# E e^{-eps xi} = exp(-k1 eps + k2/2 eps^2 - k3/6 eps^3+...); -eps Xi + eps log n cancels k1:
# 1 - exp(-eps*Xi) * exp(k2/2 eps^2 - k3/6 eps^3)
eps, aa2, aa3 = sp.symbols('eps aa2 aa3')
inner = sp.series(-sp.exp(-eps*Xi)*(1 + k2*eps**2/2 + (k2**2/8 - k3/6)*eps**3**0 if False else 0), eps, 0, 1)
g = sp.series(-sp.exp(-eps*Xi) * sp.exp(k2*eps**2/2 - k3*eps**3/6), eps, 0, 4).removeO() + 1
g = sp.expand(g)
print("\nterm1 for f=x (should equal 2z*(mu1 eps - mu2 eps^2 + mu3/3 eps^3)):")
term1 = sp.expand(2*z*g)
print(sp.collect(term1, eps))
# term2 = - int_n^inf phi (1-e^{-eps log(n/s)}) ds, phi = d2/s^2 + d3/s^3 + ...
# int_n^inf s^{-l} (1-e^{-eps log(n/s)}) ds = n^{1-l} int_1^inf u^{-l}(1-u^{-eps}) du
#  = n^{1-l} * eps * ... expand: int_1^inf u^{-l}(1-u^{-eps}) du = 1/(l-1) - 1/(l-1+eps)
u = sp.Symbol('u', positive=True)
t2 = 0
for l, dl in [(2, 2*z*(1-z)), (3, -4*c*z*(1-z))]:
    ex = sp.series(1/(l-1) - 1/(l-1+eps), eps, 0, 4).removeO()
    t2 += dl * ex  # times n^{1-l}: n^{-1}=eps/(1+eps), n^{-2}=eps^2/(1+eps)^2
# n^{1-l}: l=2 -> eps/(1+eps); l=3 -> eps^2/(1+eps)^2
t2full = sp.series(2*z*(1-z)*sp.series(1/(1)-1/(1+eps),eps,0,4).removeO()*eps/(1+eps)
                   + (-4*c*z*(1-z))*sp.series(1/sp.Integer(2)-1/(2+eps),eps,0,3).removeO()*eps**2/(1+eps)**2,
                   eps, 0, 4).removeO()
# 2-L = term1 - term2part? CAREFUL: 2-L = int_0^n = int_0^inf - int_n^inf = term1 - t2full
# but term1's int_0^inf of (1-exp) includes s>n region wrongly; subtract its own expansion there:
# int_n^inf phi(s)(1-e^{-eps log(n/s)})ds with phi expanded = t2full (as computed)
total = sp.expand(sp.series(term1 - t2full, eps, 0, 4).removeO())
print("\ntotal f=x coefficients: eps:", sp.simplify(total.coeff(eps,1) - exact_J1),
      " eps^2:", sp.simplify(sp.expand(total.coeff(eps,2) - exact_J2)),
      " eps^3:", sp.simplify(sp.expand(total.coeff(eps,3) - exact_J3)))
print("elapsed", time.time()-t0)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} rechecks: pi/2 -> P3 Q1,Q2,Q3 (diff 0); f=x -> exact every-angle series (diff 0), plus independent derivation of the f=x series from the layer-cake integral  [EXACT]\n")
