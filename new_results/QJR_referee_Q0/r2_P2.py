# Referee QJR: P2c and P2a, substitutions, antiderivatives, branch/continuation checks (own code).
from rcommon import *
x,phi,k,tau,w,beta_=sp.symbols('x phi k tau w beta',positive=True)
qs=sp.symbols('q',positive=True)
c=sp.sqrt(sp.Rational(2,3))
print('=== P2c')
B=(x**2+1)/(2*c*x)
print(' beta(1)=sqrt(3/2):',sp.simplify(B.subs(x,1)-sp.sqrt(sp.Rational(3,2)))==0,' beta(sqrt3)=sqrt2:',sp.simplify(B.subs(x,sp.sqrt(3))-r2)==0)
print(' 2beta^2-3 == 3(x^2-1)^2/(4x^2):',sp.simplify(2*B**2-3-3*(x**2-1)**2/(4*x**2))==0)
print(' dbeta/sqrt(2beta^2-3) == dx/(sqrt2 x) (x>1):',sp.simplify(sp.diff(B,x)/(sp.sqrt(3)*(x**2-1)/(2*x))-1/(r2*x))==0)
print(' (beta+1)/(beta-1) == (x^2+2cx+1)/(x^2-2cx+1):',sp.simplify((B+1)/(B-1)-(x**2+2*c*x+1)/(x**2-2*c*x+1))==0)
# P2c integrand in q -> in beta: 2*asinh(1/q)/beta * q dq/sqrt(2q^2-1) with q dq = beta dbeta
bb=sp.sqrt(1+qs**2)
print(' asinh(1/q)=1/2 log((beta+1)/(beta-1)) (num, q=0.8):',sp.N(sp.asinh(1/qs)-sp.log((bb+1)/(bb-1))/2,50).subs(qs,sp.Rational(4,5)).evalf(50))
# antiderivative: d/dx Re Li2(x w), w=e^{i alpha_c}=(sqrt2+i)/sqrt3; Re log(1-xw) = 1/2 log(x^2-2cx+1)
om=(r2+sp.I)/sp.sqrt(3)
print(' |1-x w|^2 == x^2-2cx+1:',sp.simplify(sp.expand((1-x*om)*(1-x*sp.conjugate(om)))-(x**2-2*c*x+1))==0,
      ' |1+x w|^2 == x^2+2cx+1:',sp.simplify(sp.expand((1+x*om)*(1+x*sp.conjugate(om)))-(x**2+2*c*x+1))==0)
print(' Im(x w)=x/sqrt3 != 0 for x>0, so x*w, -x*w never on [1,oo): Li2 analytic, d/dx Li2(xw) = -log(1-xw)/x')
mp.mp.dps=110
s2=mp.sqrt(2); A=1+s2; Bv=s2-1; Lv=mp.log(A); lv=mp.log(2); l3v=mp.log(3); alv=mp.atan(1/s2)
Li=lambda z: mp.polylog(2,z)
chi=lambda z:(Li(z)-Li(-z))/2
omv=mp.mpc(s2,1)/mp.sqrt(3)
G=lambda xx: 2*s2*mp.re(chi(xx*omv))
# numeric derivative check of G at interior points
for xx in [mp.mpf('1.1'),mp.mpf('1.5'),mp.mpf('1.7')]:
    integrand=(mp.log(xx**2+2*mp.sqrt(mp.mpf(2)/3)*xx+1)-mp.log(xx**2-2*mp.sqrt(mp.mpf(2)/3)*xx+1))/(s2*xx)
    print('   G\'(x)-integrand at x=%s: %s'%(mp.nstr(xx,3),mp.nstr(mp.diff(G,xx)-integrand,3)))
P2c_quad_x=mp.quad(lambda xx:(mp.log(xx**2+2*mp.sqrt(mp.mpf(2)/3)*xx+1)-mp.log(xx**2-2*mp.sqrt(mp.mpf(2)/3)*xx+1))/(s2*xx),[1,mp.sqrt(3)])
P2c_G=G(mp.sqrt(3))-G(1)
print(' P2c: x-form quad vs G(sqrt3)-G(1):',mp.nstr(abs(P2c_quad_x-P2c_G),3))
# steps (a)-(e)
z=mp.mpc(s2,1); t=1/z; y=(1-t)/(1+t)
print(' (a) Re chi2(e^{i a_c}) - (pi^2/8 - pi a_c/4):',mp.nstr(abs(mp.re(chi(omv))-(mp.pi**2/8-mp.pi*alv/4)),3))
print(' (b) t=(sqrt2-i)/3:',mp.nstr(abs(t-mp.mpc(s2,-1)/3),3),'  Re chi2(z)+Re chi2(t) - (pi^2/4 - pi a_c/2):',mp.nstr(abs(mp.re(chi(z)+chi(t))-(mp.pi**2/4-mp.pi*alv/2)),3))
print('     log z, log(-z) principal:',mp.nstr(mp.log(z),12),mp.nstr(mp.log(-z),12))
print(' (c) y=(1-t)/(1+t) vs b e^{i pi/4}:',mp.nstr(abs(y-Bv*mp.expjpi(mp.mpf(1)/4)),3),
      ' log((1+t)/(1-t)) vs L - i pi/4:',mp.nstr(abs(mp.log((1+t)/(1-t))-mp.mpc(Lv,-mp.pi/4)),3))
print('     Re chi2(t)+Re chi2(y) - (pi^2/8 - L log3/4 - pi a_c/8):',mp.nstr(abs(mp.re(chi(t)+chi(y))-(mp.pi**2/8-Lv*l3v/4-mp.pi*alv/8)),3))
print(' (d) Re chi2(y) - (Re Li2(y) - Li2(-b^4)/16):',mp.nstr(abs(mp.re(chi(y))-(mp.re(Li(y))-Li(-Bv**4)/16)),3))
print(' (e) Re Li2(b e^{i pi/4}) - (Li2(b)/2+pi^2/192):',mp.nstr(abs(mp.re(Li(y))-(Li(Bv)/2+mp.pi**2/192)),3))
# Continuation of chi2-Landen: test the COMPLEX identity on a grid over C\((-oo,0]U[1,oo)), incl. near the cuts
def landen_chi(tt): return chi(tt)+chi((1-tt)/(1+tt))-(mp.pi**2/8+mp.log(tt)*mp.log((1+tt)/(1-tt))/2)
mx=0; worst=None
mp.mp.dps=40
pts=[]
for re_ in ['-5','-1.5','-1','-0.9','-0.5','-0.01','0.01','0.3','0.5','0.99','1','1.01','2','7']:
    for im_ in ['1e-8','-1e-8','0.01','-0.01','0.5','-0.5','1','-3','10']:
        pts.append(mp.mpc(re_,im_))
pts.append(mp.mpc(s2,-1)/3)
for tt in pts:
    d=abs(landen_chi(tt))
    if d>mx: mx=d; worst=tt
print(' chi2-Landen as a complex identity on %d points of C\\((-oo,0]U[1,oo)) (incl. 1e-8 from both cuts): max |residual| = %s'%(len(pts),mp.nstr(mx,3)))
# What fails if one leaves the domain (sanity: residual on the cut sides differs is expected)
print('   on the cut t=-0.5 (principal values): residual =',mp.nstr(abs(landen_chi(mp.mpf('-0.5'))),5),'(not used)')
mp.mp.dps=110
# Combine (a)-(e) exactly as claimed
P2c_claim=s2*(Li(Bv)-Li(-Bv**4)/8+Lv*l3v/2+mp.pi**2/96-mp.pi*alv/4)
print(' P2c = 2sqrt2[Re chi2(sqrt2+i)-Re chi2(w)] vs Lemma P closed form:',mp.nstr(abs(P2c_G-P2c_claim),3))
# symbolic solve of (b)-(e)
X,Tt,Y,ReLy,Lm,Mm,PI,AL,LL,L3=sp.symbols('X T Y ReLiy Lam M PI AL LL L3')
sol=sp.solve([sp.Eq(X+Tt,PI**2/4-PI*AL/2),sp.Eq(Tt+Y,PI**2/8-LL*L3/4-PI*AL/8),sp.Eq(Y,ReLy-Mm/16),sp.Eq(ReLy,Lm/2+PI**2/192)],[X,Tt,Y,ReLy])
P2c_sym=sp.expand(2*r2*(sol[X]-(PI**2/8-PI*AL/4)))
print(' symbolic P2c from (a)-(e):',P2c_sym, '  equals printed:',sp.expand(P2c_sym-r2*(Lm-Mm/8+LL*L3/2+PI**2/96-PI*AL/4))==0)

print('=== P2a')
q_=sp.sinh(phi)
cc=c
Gp=sp.diff(sp.acosh(cc*sp.cosh(phi)),phi)
print(' c^2cosh^2-1 == (2sinh^2-1)/3:',sp.simplify((cc*sp.cosh(phi))**2-1-(2*sp.sinh(phi)**2-1)/3)==0)
for pv in ['0.7','0.8','0.88']:
    pv=sp.Float(pv,60)
    print('   G\'(phi)-sqrt2 sinh/sqrt(2sinh^2-1) at phi=%s:'%sp.N(pv,3),sp.N((Gp-r2*sp.sinh(phi)/sp.sqrt(2*sp.sinh(phi)**2-1)).subs(phi,pv),50))
phi0=sp.asinh(1/r2)
print(' G(phi0)=acosh(c cosh phi0):',sp.simplify(cc*sp.cosh(phi0)),'  c cosh L =',sp.simplify(cc*sp.sqrt(2)),' acosh(2/sqrt3)-log3/2 =',sp.N(sp.acosh(2/sp.sqrt(3))-sp.log(3)/2,60))
# W in five forms numerically
c_=mp.sqrt(mp.mpf(2)/3); ph0=mp.asinh(1/s2)
W1=mp.quad(lambda p:mp.acosh(c_*mp.cosh(p)),[ph0,Lv])
W2=mp.quad(lambda kk:mp.acosh(kk/mp.sqrt(1-kk**2))/kk,[1/s2,c_])
W3=mp.quad(lambda tt:mp.atanh(tt)*tt/(2-tt**2),[0,1/s2])
W4=mp.quad(lambda tt:mp.log(1+tt)*(1/(s2-tt)-1/(s2+tt))/4,[-1/s2,0,1/s2])
W5=mp.quad(lambda ww:mp.log(ww)*(1/(A-ww)-1/(ww+Bv))/4,[Bv/s2,1,A/s2])
print(' W forms (110 dps): |W2-W1|=%s |W3-W1|=%s |W4-W1|=%s |W5-W1|=%s'%tuple(mp.nstr(abs(v-W1),3) for v in [W2,W3,W4,W5]))
# Direct P2a by quadrature in phi, then IBP form
P2a_phi=2*mp.quad(lambda p:p*mp.sinh(p)/mp.sqrt(2*mp.sinh(p)**2-1),[ph0,Lv])
print(' P2a phi-form vs sqrt2[L log3/2 - W]:',mp.nstr(abs(P2a_phi-s2*(Lv*l3v/2-W1)),3))
# Tonelli step symbolic pieces
mu=sp.sqrt(1-k**2)/k
print(' k^2cosh^2-1 == k^2(sinh^2-mu^2):',sp.simplify(k**2*sp.cosh(phi)**2-1-k**2*(sp.sinh(phi)**2-mu**2))==0)
inner=sp.acosh(sp.sinh(phi)/mu)/k
print(' d/dphi [acosh(sinh/mu)/k] == cosh/sqrt(k^2cosh^2-1) (num at k=0.75,phi=0.8):',
      sp.N((sp.diff(inner,phi)-sp.cosh(phi)/sp.sqrt(k**2*sp.cosh(phi)**2-1)).subs({k:sp.Rational(3,4),phi:sp.Rational(4,5)}),50))
kt=(2-tau**2)**sp.Rational(-1,2)
print(' tau-substitution: (k/sqrt(1-k^2))^2=1/(1-tau^2) (both sides >0 for 0<=tau<1):',sp.simplify((kt**2/(1-kt**2))-1/(1-tau**2))==0,
      ' dk/k = tau dtau/(2-tau^2):',sp.simplify(sp.diff(kt,tau)/kt-tau/(2-tau**2))==0,
      ' tau/(2-tau^2)=K/2:',sp.simplify(tau/(2-tau**2)-(1/(r2-tau)-1/(r2+tau))/2)==0)
# F_- , F_+ derivatives
Fm=-sp.log(w)*sp.log(1-w/a)-sp.polylog(2,w/a)
Fp=sp.log(w)*sp.log(1+w/b)+sp.polylog(2,-w/b)
dFm=sp.diff(Fm,w).rewrite(sp.log); dFp=sp.diff(Fp,w).rewrite(sp.log)
print(' F_-\' raw:',dFm); print(' F_+\' raw:',dFp)
# polylog'(2,z) = polylog(1,z)/z = -log(1-z)/z ; check algebraically after this replacement
dFm2=sp.diff(Fm,w).replace(sp.polylog,lambda s_,z_: -sp.log(1-z_) if s_==1 else sp.polylog(s_,z_))
dFp2=sp.diff(Fp,w).replace(sp.polylog,lambda s_,z_: -sp.log(1-z_) if s_==1 else sp.polylog(s_,z_))
e1=sp.expand_log(sp.simplify(dFm2-sp.log(w)/(a-w)),force=True); e2=sp.expand_log(sp.simplify(dFp2-sp.log(w)/(w+b)),force=True)
print(' F_-\' - log w/(a-w) =',sp.simplify(e1),'   F_+\' - log w/(w+b) =',sp.simplify(e2))
for wv in ['0.4','1.0','1.6']:
    print('   numeric at w=%s:'%wv,sp.N((dFm2-sp.log(w)/(a-w)).subs(w,sp.Rational(wv)),40),sp.N((dFp2-sp.log(w)/(w+b)).subs(w,sp.Rational(wv)),40))
# exact evaluation of W with logS and Li2 symbols
def Fm_e(wv): return -logS(wv)*logS(1-wv/a)-Li2(wv/a)
def Fp_e(wv): return logS(wv)*logS(1+wv/b)+Li2(-wv/b)
lo=b/r2; hi=a/r2
print(' arguments: w/a at ends:',canon(lo/a),canon(hi/a),' -w/b:',canon(-lo/b),canon(-hi/b))
Wex=sp.expand(((Fm_e(hi)-Fm_e(lo))-(Fp_e(hi)-Fp_e(lo)))/4)
Wpr=lam*L/2-L*l3/2+(Li2(-1/r2)-Li2(1/r2)+Li2(b**2/r2)-Li2(-a**2/r2))/4
print(' W exact =',Wex)
print(' W exact - printed W = ',sp.simplify(Wex-Wpr))
print(' |W exact (num) - W1| =',mp.nstr(abs(numval(Wex)-W1),3))
import json
tv={kk:sp.sympify(vv) for kk,vv in json.load(open('r3_table_values.json')).items()}
D15=M-10*Lam-L**2-lam*L+7*pi**2/8
br=sp.expand(tv[str(canon(-1/r2))]-tv[str(canon(1/r2))]+D15)
print(' bracket via table (2),(3),(15):',br,'  equals printed -12Lam+M-2L^2+5pi^2/6:',sp.expand(br-(-12*Lam+M-2*L**2+sp.Rational(5,6)*pi**2))==0)
Wred=sp.expand(lam*L/2-L*l3/2+br/4)
P2a_red=sp.expand(r2*(L*l3/2-Wred))
P2a_pr=r2*(3*Lam-M/4-lam*L/2+L*l3+L**2/2-5*pi**2/24)
print(' P2a reduced:',P2a_red,'  equals Lemma P:',sp.expand(P2a_red-P2a_pr)==0)
