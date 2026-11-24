# Referee QJR: direct quadrature of the ORIGINAL definitions (own code), compared with closed forms.
import mpmath as mp, sys, json, time
DPS=int(sys.argv[1]) if len(sys.argv)>1 else 50
mp.mp.dps=DPS+15
t0=time.time()
s2=mp.sqrt(2); a=1+s2; b=s2-1; L=mp.log(a); lam=mp.log(2); l3=mp.log(3); al=mp.atan(1/s2); pi=mp.pi
Li=lambda x: mp.polylog(2,x)
CF={'P1': 3*s2*Li(b)-lam**2/2-s2/2*lam*L+(1+3*s2/2)*L**2-(mp.mpf(1)/48+3*s2/16)*pi**2,
    'P2a': s2*(3*Li(b)-Li(-b**4)/4-lam*L/2+L*l3+L**2/2-5*pi**2/24),
    'P2c': s2*(Li(b)-Li(-b**4)/8+L*l3/2+pi**2/96-pi*al/4)}
CF['Q0']=-pi**2/8+s2/16*pi**2+3*s2/2*pi*al-3*lam**2-9*s2*L*l3+(6+6*s2)*L**2-6*s2*Li(b)+9*s2/4*Li(-b**4)
def q_of(s): return mp.sqrt(s**2+(1-s)**2)
def f1(s):   # P1 integrand, original definition; 1-q^2 = 2s(1-s) exactly
    with mp.workdps(2*mp.mp.dps):
        q=q_of(s); return (2*mp.log(1+q)-(1+q)*mp.log(2*q))/(2*s*(1-s))
def f2a(s): q=q_of(s); return mp.asinh(q)/mp.sqrt(1+q**2)
def f2c(s): q=q_of(s); return mp.asinh(1/q)/mp.sqrt(1+q**2)
res={}
which=sys.argv[2] if len(sys.argv)>2 else 'P'
if 'P' in which:
    for name,f in [('P1',f1),('P2a',f2a),('P2c',f2c)]:
        v=2*mp.quad(f,[0,mp.mpf(1)/8,mp.mpf(1)/4,mp.mpf(1)/2])   # q(s)=q(1-s)
        res[name]=v
        print(name,'quad =',mp.nstr(v,DPS),' |quad-closed| =',mp.nstr(abs(v-CF[name]),3),flush=True)
    v=6*(res['P1']-res['P2a']-res['P2c'])
    print('6(P1-P2a-P2c) quad - Q0 closed =',mp.nstr(abs(v-CF['Q0']),3),flush=True)
if 'J' in which:   # Q0 = 6 int_0^1 J(q(s)) ds with J by quadrature of its definition
    def J(q):
        g=lambda r:(1/(1+r+mp.sqrt(q**2+r**2))-1/(1+q))/(r*(1+r))
        return mp.quad(g,[0,q/4,q,1,4,40,mp.inf])
    v=12*mp.quad(lambda s:J(q_of(s)),[0,mp.mpf(1)/8,mp.mpf(1)/4,mp.mpf(1)/2])
    print('Q0 via 6 int J (nested quad) =',mp.nstr(v,DPS),' |diff| =',mp.nstr(abs(v-CF['Q0']),3),flush=True)
if 'K' in which:   # Q0 = -pi^{-1/2} int_0^oo k0(y)^3 dy
    def k0(y):
        g=lambda x:-mp.expm1(-(x/y)**2)*mp.exp(-x)/x
        pts=[0,y,1,mp.inf] if y<1 else [0,1,y,mp.inf]
        return mp.quad(g,pts)
    v=-mp.quad(lambda y:k0(y)**3,[0,mp.mpf(1)/100,mp.mpf(1)/10,1,10,mp.inf])/mp.sqrt(pi)
    print('Q0 via -pi^-1/2 int k0^3 =',mp.nstr(v,DPS),' |diff| =',mp.nstr(abs(v-CF['Q0']),3),flush=True)
print('Q0 closed form =',mp.nstr(CF['Q0'],DPS),'  time %.0fs'%(time.time()-t0))
