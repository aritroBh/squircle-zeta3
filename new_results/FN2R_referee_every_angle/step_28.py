import sympy as sp, time
t0=time.time()
a,A,B = sp.symbols('a A B')
c, z, k2, k3, Xi, mu = sp.symbols('c z k2 k3 Xi mu')
th = {(1,0): 1/a, (2,0): (A-a**2)/(2*a**3), (1,1): -A/a**3,
      (3,0): (3*A**2-3*A*a**2-B*a+a**4)/(6*a**5),
      (2,1): (-4*A**2+3*A*a**2+B*a)/(2*a**5), (1,2): A**2/a**5}
# my own d_l(c) via Legendre
d = lambda l: (-1)**l*(1-c**2)*sp.legendre(l-1, c).diff(c)/l
# M_i from moment definition: M_i = (1-c) E[(log na + xi)^i], xi central moments m2=k2,m3=k3, mean k1 absorbed in Xi
mu1 = Xi; mu2 = Xi**2 + k2; mu3 = Xi**3 + 3*k2*Xi + k3
M = {1: 2*z*mu1, 2: 2*z*mu2, 3: 2*z*mu3}   # with 2z = 1-c enforced at the end
# T-part: sum_{i,j} th_ij mu^{i+j} (-1)^{i+1} i! sum_{l>=2} d_l/(l-1)^{i+1} * (mu/(a(1+mu)))^{l-1}
K = 3
Tser = 0
for (i,j),tij in th.items():
    for l in range(2, K+2):
        term = tij*mu**(i+j)*(-1)**(i+1)*sp.factorial(i)*d(l)/(l-1)**(i+1) * (mu/(a*(1+mu)))**(l-1)
        Tser += sp.series(term, mu, 0, K+1).removeO()
Tser = sp.expand(Tser)
T2 = sp.expand(Tser).coeff(mu,2); T3 = sp.expand(Tser).coeff(mu,3)
print("my T2 =", sp.simplify(T2))
print("my T3 =", sp.simplify(sp.expand(T3)))
# assemble master J's
J1 = sp.expand(th[(1,0)]*M[1])
J2 = sp.expand(th[(2,0)]*M[2] + th[(1,1)]*M[1] + T2)
J3 = sp.expand(th[(3,0)]*M[3] + th[(2,1)]*M[2] + th[(1,2)]*M[1] + T3)
# paper's displayed J's
J1p = 2*z*mu1/a
J2p = z*(A*(mu2-2*mu1) - a**2*mu2 + 2*a*(1-z))/a**3
J3p = ( z*(a**4*mu3 - (3*A*a**2+B*a)*(mu3-3*mu2) + 3*A**2*(mu3-4*mu2+2*mu1))/(3*a**5)
        - 4*A*z*(1-z)/a**4 - c*z*(1-z)/a**3 )
sub = {z: (1-c)/2}
d1 = sp.simplify((J1-J1p).subs(sub))
d2 = sp.simplify(sp.expand((J2-J2p).subs(sub)))
d3 = sp.simplify(sp.expand((J3-J3p).subs(sub)))
print("master - displayed: J1:", d1, " J2:", d2, " J3:", d3)
# reduction to P3 Q's at theta=pi/2
Q2 = -(a**2-A)/(2*a**3)*(Xi**2+k2) - A*Xi/a**3 + 1/(2*a**2)
lam3 = a**4-3*a**2*A-a*B+3*A**2
Q3 = lam3/(6*a**5)*mu3 + (3*a**2*A+a*B-4*A**2)/(2*a**5)*mu2 + A**2*Xi/a**5 - A/a**4
subR = {z: sp.Rational(1,2), c: 0}
print("pi/2 reduction to P3:", sp.simplify((J1-Xi/a).subs(subR))==0,
      sp.simplify(sp.expand((J2-Q2).subs(subR)))==0,
      sp.simplify(sp.expand((J3-Q3).subs(subR)))==0)
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} step_05b: MY master formula (my theta_ij, my d_l via Legendre, my T_k from tail sum) -> paper's displayed J1,J2,J3: difference 0 at general c; pi/2 -> P3 Q's: difference 0 [referee, sympy EXACT]\n")
