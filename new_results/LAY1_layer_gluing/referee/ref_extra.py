# Extra checks: Lemma 5 at (1.5,0.6) with a longer tail; h+h'' min at tau=0.01; R3/tau^4 at 1/tau=10.
import sys
from mpmath import mp, mpf, quad, exp, gamma, inf, nstr, cos, sin, diff, pi, zeta, sqrt
mp.dps = 60
lam, mu = mpf('1.5'), mpf('0.6')
G = lambda l,m: quad(lambda s: exp(l*s)*((1+exp(-2*s))**m-1), [0,1,2,4,8,16,32,64,128,256,512,inf])
R = (gamma(-lam/2)*gamma(lam/2-mu)/gamma(-mu)+2/lam+2/(2*mu-lam))/2
print("Lemma 5 at (1.5,0.6), tail to 512: |LHS-RHS| =", nstr(abs(G(lam,mu)+G(2*mu-lam,mu)-R),3))
sys.argv=['x','60']; exec(open('ref5_endtoend.py').read().split('z = zeta')[0])
tau = mpf(1)/10
I,_ = I_squircle(tau); print("squircle 1/tau=10 R3/tau^4 =", nstr((I-4*sqrt(2)-sqrt(2)*zeta(2)*tau**2-sqrt(2)*zeta(3)*tau**3)/tau**4,20))
I,_ = I_softmax(tau);  print("softmax  1/tau=10 R3/tau^4 =", nstr((I-8-pi**2/3*tau**2-2*zeta(3)*tau**3)/tau**4,20))
mp.dps = 300
V=[(1,1),(1,-1),(-1,1),(-1,-1)]
def h(th,tau):
    d=[v[0]*cos(th)+v[1]*sin(th) for v in V]; m=max(d); return m*sum(exp((x-m)/tau) for x in d)**tau
for tau in (mpf('0.01'),):
    k = lambda th: h(th,tau)+diff(lambda t: h(t,tau), th, 2)
    print("tau=0.01: h+h'' at pi/4 =", nstr(k(pi/4),8), " at pi/4-0.05 =", nstr(k(pi/4-mpf('0.05')),8), " at 0.3 =", nstr(k(mpf('0.3')),8))
    mp.dps = 55
    print("   same at dps 55 (pi/4):", nstr(h(pi/4,tau)+diff(lambda t: h(t,tau), pi/4, 2),8))
