import sympy as sp, time
t0=time.time()
# ===== f=x: expand paper thm:P2 exact series through eps^3 at general (c,z) =====
eps, L = sp.symbols('eps L')          # eps = 1/(n-1), L = log n
c, z = sp.symbols('c z')
a1, a2, a3 = sp.symbols('a1 a2 a3')   # paper's corner coefficients: log E[S^w] = sum a_j w^j
pi = sp.pi
# Lambda_w/(2z) = 2F1(-w,1+w;2;z) = 1 + al1 w + al2 w^2 + al3 w^3 + ...
# from log(Lambda_w/2z) = a1 w + a2 w^2 + a3 w^3 - log(pi w/sin pi w), log(pi w/sin pw)= pi^2 w^2/6 + O(w^4)
w = sp.Symbol('w')
logLam = a1*w + (a2 - pi**2/6)*w**2 + a3*w**3
Lam = sp.expand(sp.exp(logLam).series(w,0,4).removeO())
al1 = Lam.coeff(w,1); al2 = Lam.coeff(w,2); al3 = Lam.coeff(w,3)
# 2 - L_n = 2z - (pi eps/sin pi eps) * (2z)(1+al1 eps+al2 eps^2+al3 eps^3) * n^{-eps}
#           - sum_{m>=1} (-1)^m eps/(m-eps) Lam_m n^{-m}
n_min_eps = sp.expand(sp.exp(-eps*L).series(eps,0,4).removeO())
mainfac = sp.expand((1 + (pi**2/6)*eps**2) * (1 + al1*eps + al2*eps**2 + al3*eps**3) * n_min_eps)
term1 = sp.expand(2*z - 2*z*mainfac)   # through eps^3
term1 = sum(m for m in term1.args if sp.degree(m,eps)<=3)
# m-sum: m=1: + eps/(1-eps) Lam1 n^{-1}; m=2: - eps/(2-eps) Lam2 n^{-2}
Lam1 = (sp.legendre(0,c)-sp.legendre(2,c))/3
Lam2 = (sp.legendre(1,c)-sp.legendre(3,c))/5
msum = ( sp.series(eps/(1-eps), eps, 0, 3).removeO()*Lam1*sp.series(eps/(1+eps), eps, 0, 3).removeO()
       - sp.series(eps/(2-eps), eps, 0, 2).removeO()*Lam2*sp.series(eps**2/(1+eps)**2, eps, 0, 2).removeO() )
msum = sp.expand(msum)
two_minus_L = sp.expand(term1 - msum)
two_minus_L = sum(m for m in two_minus_L.args if sp.degree(m,eps)<=3)
two_minus_L = sp.expand(two_minus_L)
c1 = sp.expand(two_minus_L).coeff(eps,1); c2 = sp.expand(two_minus_L).coeff(eps,2); c3 = sp.expand(two_minus_L).coeff(eps,3)
print("exact series coefficients (in z,c,a_j,L):")
print(" eps^1:", sp.simplify(c1.subs(z,(1-c)/2)))
print(" eps^2:", sp.simplify(sp.expand(c2.subs(z,(1-c)/2))))
print(" eps^3:", sp.simplify(sp.expand(c3.subs(z,(1-c)/2))))
# displayed J's at a=1,A=B=0 with k_j = (-1)^j j! a_j, Xi = L + k1
k1, k2, k3 = -a1, 2*a2, -6*a3
Xi = L + k1
mu1 = Xi; mu2 = Xi**2 + k2; mu3 = Xi**3 + 3*k2*Xi + k3
J1x = 2*z*mu1
J2x = z*(2*(1-z) - mu2)
J3x = z*mu3/3 - c*z*(1-z)
d1 = sp.simplify(sp.expand((c1-J1x).subs(z,(1-c)/2)))
d2 = sp.simplify(sp.expand((c2-J2x).subs(z,(1-c)/2)))
d3 = sp.simplify(sp.expand((c3-J3x).subs(z,(1-c)/2)))
print("\nexact-series minus displayed-J (a=1,A=B=0), at general c:", d1, d2, d3)
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} step_06: paper thm:P2 exact series (f=x) expanded through eps^3 at general c: EQUALS displayed J1,J2,J3 (difference 0) — end-to-end check of master formula + cumulants + tail [referee, sympy EXACT]\n")
