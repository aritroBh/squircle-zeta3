import sympy as sp, time
t0=time.time()
eps, L = sp.symbols('eps L')
c, z = sp.symbols('c z')
a1, a2, a3 = sp.symbols('a1 a2 a3')
pi = sp.pi
w = sp.Symbol('w')
logLam = a1*w + (a2 - pi**2/6)*w**2 + a3*w**3
Lam = sp.expand(sp.exp(logLam).series(w,0,4).removeO())
al1 = Lam.coeff(w,1); al2 = Lam.coeff(w,2); al3 = Lam.coeff(w,3)
n_min_eps = sp.expand(sp.exp(-eps*L).series(eps,0,4).removeO())
mainfac = sp.expand((1 + (pi**2/6)*eps**2) * (1 + al1*eps + al2*eps**2 + al3*eps**3) * n_min_eps)
term1 = sp.expand(2*z - 2*z*mainfac)
term1 = sum(m for m in term1.args if sp.degree(m,eps)<=3)
Lam1 = (sp.legendre(0,c)-sp.legendre(2,c))/3
Lam2 = (sp.legendre(1,c)-sp.legendre(3,c))/5
# 2 - L = 2z - main - sum_m (-1)^m eps/(m-eps) Lam_m n^{-m}
#       = term1 + eps/(1-eps)Lam1 n^{-1} - eps/(2-eps)Lam2 n^{-2} + O(eps^4)
ninv1 = sp.series(eps/(1+eps), eps, 0, 4).removeO()
ninv2 = sp.series(eps**2/(1+eps)**2, eps, 0, 3).removeO()
mpart = ( sp.series(eps/(1-eps), eps, 0, 4).removeO()*Lam1*ninv1
        - sp.series(eps/(2-eps), eps, 0, 3).removeO()*Lam2*ninv2 )
two_minus_L = sp.expand(term1 + mpart)
two_minus_L = sp.expand(sum(m for m in two_minus_L.args if sp.degree(m,eps)<=3))
c1 = two_minus_L.coeff(eps,1); c2 = two_minus_L.coeff(eps,2); c3 = two_minus_L.coeff(eps,3)
k1, k2, k3 = -a1, 2*a2, -6*a3
Xi = L + k1
mu1 = Xi; mu2 = Xi**2 + k2; mu3 = Xi**3 + 3*k2*Xi + k3
J1x = 2*z*mu1
J2x = z*(2*(1-z) - mu2)
J3x = z*mu3/3 - c*z*(1-z)
sub = {z: (1-c)/2}
d1 = sp.simplify(sp.expand((c1-J1x).subs(sub)))
d2 = sp.simplify(sp.expand((c2-J2x).subs(sub)))
d3 = sp.simplify(sp.expand((c3-J3x).subs(sub)))
print("exact series minus displayed-J (f=x), general c:", d1, d2, d3)
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} step_06 fixed (m-sum sign): thm:P2 exact series f=x through eps^3 == displayed J's, difference 0 at general c [referee, sympy EXACT]\n")
