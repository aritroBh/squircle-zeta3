# Referee QJR: Q0 = -pi^{-1/2} int_0^oo k0(y)^3 dy, k0(y)=int_0^oo (1-e^{-u^2}) e^{-y u} du/u  (x=y u), finer splitting.
import mpmath as mp, sys, time
DPS=int(sys.argv[1]); mp.mp.dps=DPS+20; t0=time.time()
s2=mp.sqrt(2); a=1+s2; b=s2-1; L=mp.log(a); lam=mp.log(2); l3=mp.log(3); al=mp.atan(1/s2); pi=mp.pi
Q0cf=-pi**2/8+s2/16*pi**2+3*s2/2*pi*al-3*lam**2-9*s2*L*l3+(6+6*s2)*L**2-6*s2*mp.polylog(2,b)+9*s2/4*mp.polylog(2,-b**4)
def k0(y):
    g=lambda u:-mp.expm1(-u*u)*mp.exp(-y*u)/u
    pts=[0,mp.mpf(1)/2,1,2,4,8]
    if y<1: pts+= [p/y for p in [mp.mpf(1)/2,1,4,16,64] if p/y>8]
    pts+= [mp.inf]
    return mp.quad(g,sorted(set(pts)) if False else pts,maxdegree=10)
ys=[0]+[mp.mpf(10)**e for e in range(-8,3)]+[mp.inf]
v=-mp.quad(lambda y:k0(y)**3,ys,maxdegree=10)/mp.sqrt(pi)
print('Q0 via -pi^-1/2 int k0^3 (u-form) =',mp.nstr(v,DPS),'\n  closed form =',mp.nstr(Q0cf,DPS),'\n  |diff| =',mp.nstr(abs(v-Q0cf),3),' time %.0fs'%(time.time()-t0))
