# Diagnose the 1.5e-17 residual of Lemma 5 at (lam,mu)=(1.5,0.6): recompute G via u-substitution K(a,mu).
from mpmath import mp, mpf, quad, exp, gamma, inf, nstr
mp.dps = 60
def G_s(l,m): return quad(lambda s: exp(l*s)*((1+exp(-2*s))**m-1), [0,1,2,4,8,16,32,64,128,256,512,inf])
def G_u(l,m):  # G = 1/2 K(-l/2, m),  K(a,m)=int_0^1 u^{a-1}((1+u)^m-1) du
    a = -l/2
    return quad(lambda u: u**(a-1)*((1+u)**m-1), [0, mpf(10)**-40, mpf(10)**-20, mpf(10)**-10, mpf(10)**-5, mpf('0.01'), 1])/2
lam, mu = mpf('1.5'), mpf('0.6')
R = (gamma(-lam/2)*gamma(lam/2-mu)/gamma(-mu)+2/lam+2/(2*mu-lam))/2
for name,Gf in (("s-integral",G_s),("u-integral",G_u)):
    a1, a2 = Gf(lam,mu), Gf(2*mu-lam,mu)
    print(name, " G(1.5,.6)=",nstr(a1,30)," G(-.3,.6)=",nstr(a2,30)," |LHS-RHS|=",nstr(abs(a1+a2-R),3))
# Closed form: K(a,mu) = (2F1(-mu, a; a+1; -1) - 1)/a  (termwise integration of the binomial series)
from mpmath import hyp2f1, mpc
def G_h(l,m):
    a = -l/2
    return (hyp2f1(-m, a, a+1, -1) - 1)/a/2
for lam,mu in ((mpf('1.5'),mpf('0.6')), (mpf('0.5'),mpf('0.2')), (mpf('-3'),mpf('-1.2')), (mpc('0.25','0.1'),mpc('0.4','-0.2'))):
    R = (gamma(-lam/2)*gamma(lam/2-mu)/gamma(-mu)+2/lam+2/(2*mu-lam))/2
    print("hyp2f1 route", nstr(lam,4), nstr(mu,4), " |LHS-RHS| =", nstr(abs(G_h(lam,mu)+G_h(2*mu-lam,mu)-R),3),
          "  |G_h - G_s| =", nstr(abs(G_h(lam,mu)-G_s(lam,mu)),3))
