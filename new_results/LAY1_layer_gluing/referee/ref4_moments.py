# Referee item 4: Lemma 5 numerics + moments by direct quadrature (own code).
from mpmath import mp, mpf, mpc, quad, exp, log, gamma, zeta, inf, sqrt, fac, nstr
mp.dps = 60
def G(lam, mu):
    f = lambda s: exp(lam*s)*((1+exp(-2*s))**mu - 1)
    return quad(f, [0,1,2,4,8,16,32,64,inf])
def rhs(lam, mu):
    return (gamma(-lam/2)*gamma(lam/2-mu)/gamma(-mu) + 2/lam + 2/(2*mu-lam))/2
pts = [(mpf('0.5'), mpf('0.2')), (mpf('-0.7'), mpf('-0.3')), (mpc('0.25','0.1'), mpc('0.4','-0.2')),
       (mpf('1.5'), mpf('0.6'))]   # last one: outside |lam|<2,|mu-lam/2|<1? |0.6-0.75|<1 inside; also try a far point
pts.append((mpf('-3'), mpf('-1.2')))  # outside the stated domain but where both G's converge
for lam,mu in pts:
    L = G(lam,mu)+G(2*mu-lam,mu); R = rhs(lam,mu)
    print("lam,mu=",nstr(lam,5),nstr(mu,5)," LHS=",nstr(L,20)," |LHS-RHS|=",nstr(abs(L-R),3))
g = lambda s: log(1+exp(-2*s))
Q = lambda f: quad(f,[0,0.5,1,2,4,8,16,32,64,inf])
g01 = Q(g); g02 = Q(lambda s: g(s)**2)/2; g11 = Q(lambda s: s*g(s))
g03 = Q(lambda s: g(s)**3)/6; g12 = Q(lambda s: s*g(s)**2)/2; g21 = Q(lambda s: s**2*g(s))/2
z = zeta
print("gamma01 - zeta(2)/4         =", nstr(g01-z(2)/4,3))
print("gamma02+gamma11 - zeta(3)/4 =", nstr(g02+g11-z(3)/4,3))
print("gamma21 - 7zeta(4)/64       =", nstr(g21-7*z(4)/64,3))
print("J=gamma12+gamma03 - zeta(4)/32 =", nstr(g12+g03-z(4)/32,3))
print("2g03+2g12+4g21 - zeta(4)/2  =", nstr(2*g03+2*g12+4*g21-z(4)/2,3))
I3 = Q(lambda s: s**3/(1+exp(2*s)))
print("int s^3/(1+e^{2s}) - 21zeta(4)/64 =", nstr(I3-21*z(4)/64,3))
# c4 from the c4 formula: per kink Phi_3 = H[g3 + g^3/6 + |s| g^2/2 - s^2 g/2] (T=H,S=0,a=0,g2=0), 4 kinks, even
for name,H,g3c in (("squircle",1/sqrt(2),mpf(-4)/3),("softmax",mpf(1),mpf(1)/3)):
    per_half = Q(lambda s: g3c*s**3/(1+exp(2*s)) + g(s)**3/6 + s*g(s)**2/2 - s**2*g(s)/2)
    c4 = 4*2*H*per_half
    claim = -33*sqrt(2)/16*z(4) if name=="squircle" else z(4)/4
    print(name,"c4 (quadrature of Phi_3) =",nstr(c4,30)," claim=",nstr(claim,30)," diff=",nstr(c4-claim,3))
# diagonal closed form
for tt in (mpf('0.3'), mpf('-0.4')):
    print("diag tau=",tt," G(t,t)-closed =", nstr(G(tt,tt)-(gamma(-tt/2)**2/gamma(-tt)+4/tt)/4,3))
