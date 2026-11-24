
import sympy as sp, pickle, mpmath as mp
WD = '<project>/research/08_theory/K1f_hexagon_order5'
mp.mp.dps = 50

# ---- test edge data ----
Bq=sp.Rational(1,3); rq=sp.Rational(2); kapq=rq-1
A1q=1-Bq; A0q=A1q-Bq*kapq
L2q=sp.Rational(5); q0q=sp.Rational(1); V02q=sp.Rational(3)
lA0q=sp.Rational(31,100); lA1q=-sp.Rational(17,100); lBq=sp.Rational(42,100)
print("A0,A1,B,kap =",A0q,A1q,Bq,kapq,"  A0=1-Br:",1-Bq*rq)

# ---- symbolic E5 evaluated at t=0.3, X=log t, Y=log(1-t), W=log(1+kap t) ----
A0,A1,B,kap=sp.symbols('A0 A1 B kap'); L2,q0,V02=sp.symbols('L2 q0 V02')
lA0,lA1,lB=sp.symbols('lA0 lA1 lB'); t,X,Y,W=sp.symbols('t X Y W')
local={'A0':A0,'A1':A1,'B':B,'kap':kap,'L2':L2,'q0':q0,'V02':V02,'lA0':lA0,'lA1':lA1,'lB':lB,'t':t,'X':X,'Y':Y,'W':W}
Ea=pickle.load(open(f'{WD}/E5_general_approachA.pkl','rb'))['E']
tq=sp.Rational(3,10)
subs={A0:A0q,A1:A1q,B:Bq,kap:kapq,L2:L2q,q0:q0q,V02:V02q,lA0:lA0q,lA1:lA1q,lB:lBq,
      t:tq,X:sp.log(tq),Y:sp.log(1-tq),W:sp.log(1+kapq*tq)}
E5val=sp.N(sp.sympify(Ea[5],locals=local).subs(subs),45)
print("symbolic E5(0.3) =", E5val)

# ---- independent mpmath jet d-expansion ----
NJ=6
def jmul(a,b): return [mp.fsum(a[i]*b[k-i] for i in range(k+1)) for k in range(NJ)]
def jexp(a):
    x=[mp.mpf(0)]+list(a[1:]); res=[mp.mpf(1)]+[mp.mpf(0)]*(NJ-1); pw=list(res)
    for j in range(1,NJ):
        pw=jmul(pw,x); res=[res[k]+pw[k]/mp.factorial(j) for k in range(NJ)]
    e0=mp.exp(a[0]); return [e0*r for r in res]
def jlog(a):
    x=[mp.mpf(0)]+[a[k]/a[0] for k in range(1,NJ)]
    res=[mp.mpf(0)]*NJ; pw=[mp.mpf(1)]+[mp.mpf(0)]*(NJ-1)
    for j in range(1,NJ):
        pw=jmul(pw,x); res=[res[k]+mp.mpf((-1)**(j+1))/j*pw[k] for k in range(NJ)]
    res[0]=mp.log(a[0]); return res
def jinv(a):
    r=[1/a[0]]+[mp.mpf(0)]*(NJ-1)
    for k in range(1,NJ): r[k]=-mp.fsum(a[i]*r[k-i] for i in range(1,k+1))/a[0]
    return r
def jsqrt(a): return jexp([v/2 for v in jlog(a)])
def num(x): return mp.mpf(sp.N(x,50))
tt=num(tq); A0n,A1n,Bn,kn=num(A0q),num(A1q),num(Bq),num(kapq)
L2n,q0n,V02n=num(L2q),num(q0q),num(V02q); lA0n,lA1n,lBn=num(lA0q),num(lA1q),num(lBq)
a0=A0n*tt; a1=A1n*(1-tt); a2=Bn*(1+kn*tt)
LA0=mp.log(tt)+lA0n; LA1=mp.log(1-tt)+lA1n; LA2=mp.log(1+kn*tt)+lBn
def mk(a,L): return jexp([mp.log(a),L]+[mp.mpf(0)]*(NJ-2))  # a*e^{d L}
Fj=[mp.mpf(0)]*NJ
for a,L in [(a0,LA0),(a1,LA1),(a2,LA2)]:
    Fj=[x+y for x,y in zip(Fj,mk(a,L))]
dFj=[mp.mpf(0)]*NJ
for ap,L in [(A0n,LA0),(-A1n,LA1),(Bn*kn,LA2)]:
    dFj=[x+y for x,y in zip(dFj,jexp([0,L]+[mp.mpf(0)]*(NJ-2)))]
    # note: ap*e^{dL} with ap constant -> jexp([log|ap|,...]) signed; handle sign
# rebuild dFj with signs properly
dFj=[mp.mpf(0)]*NJ
for ap,L in [(A0n,LA0),(-A1n,LA1),(Bn*kn,LA2)]:
    j=jexp([0,L]+[mp.mpf(0)]*(NJ-2)); dFj=[x+ap*y for x,y in zip(dFj,j)]
Phi=jmul(dFj,jinv(Fj))
oneoverp=[mp.mpf((-1)**k) for k in range(NJ)]
Fpow=jexp([ -v for v in jmul(oneoverp,jlog(Fj))])
yyp=q0n+L2n*tt; yy=V02n+2*q0n*tt+L2n*tt**2
Phi2=jmul(Phi,Phi)
Mmet=[ (L2n if k==0 else 0) -2*yyp*Phi[k]+yy*Phi2[k] for k in range(NJ)]
integrand=jmul(Fpow, jsqrt(Mmet))
print("jet E5(0.3)      =", mp.nstr(integrand[5],45))
print("abs diff =", mp.nstr(abs(integrand[5]-mp.mpf(str(E5val))),5))
print("rel diff =", mp.nstr(abs((integrand[5]-mp.mpf(str(E5val)))/integrand[5]),5))
