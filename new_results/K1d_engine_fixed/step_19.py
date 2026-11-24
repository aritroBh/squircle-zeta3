step_41 = r'''
import pickle, sympy as sp, mpmath as mp
mp.mp.dps = 70

# ---- (1) true series engine (K1c step_34 verbatim, Phi=Fp/F) ----
def series_engine(tt, a, ap, LA, y, yp, n=3):
    def pmul(P,Q):
        R=[mp.mpf(0)]*(n+1)
        for i in range(n+1):
            for j in range(n+1-i): R[i+j]+=P[i]*Q[j]
        return R
    def pdiv(P,Q):
        R=[mp.mpf(0)]*(n+1)
        for i in range(n+1):
            R[i]=(P[i]-sum(R[j]*Q[i-j] for j in range(i)))/Q[0]
        return R
    def plog(P):
        Pp=[P[i]*i for i in range(1,n+1)]+[mp.mpf(0)]
        Q=pdiv(Pp,P); L=[mp.mpf(0)]
        for i in range(n): L.append(Q[i]/(i+1))
        return L
    def pexp(P):
        E=[mp.mpf(1)]+[mp.mpf(0)]*n
        for k in range(1,n+1): E[k]=sum(P[j]*j*E[k-j] for j in range(1,k+1))/k
        return E
    F=[mp.mpf(0)]*(n+1); Fp=[mp.mpf(0)]*(n+1)
    for aj,LAj in zip(a,LA):
        pk=mp.mpf(1)
        for k in range(n+1):
            if k: pk=pk*LAj/k
            F[k]+=aj*pk
    for aj,LAj in zip(ap,LA):
        pk=mp.mpf(1)
        for k in range(n+1):
            if k: pk=pk*LAj/k
            Fp[k]+=aj*pk
    Phi=pdiv(Fp,F)
    lg=plog(F); inv=[mp.mpf(1),mp.mpf(-1),mp.mpf(1),mp.mpf(-1)][:n+1]
    fac=pexp(pmul([mp.mpf(-2)*x for x in lg],inv))
    yp2=yp[0]**2+yp[1]**2; yyp=y[0]*yp[0]+y[1]*yp[1]; yy=y[0]**2+y[1]**2
    Q=[mp.mpf(0)]*(n+1); Q[0]=yp2
    for k in range(1,n+1): Q[k]=-2*yyp*Phi[k]
    Phisq=pmul(Phi,Phi)
    for k in range(n+1): Q[k]+=yy*Phisq[k]
    G=pmul(fac,Q)
    sq=[mp.mpf(0)]*(n+1); sq[0]=mp.sqrt(G[0])
    for k in range(1,n+1): sq[k]=(G[k]-sum(sq[j]*sq[k-j] for j in range(1,k)))/(2*sq[0])
    return sq[3]

tt=mp.mpf('0.3')
a=[tt/2,3*(1-tt)/4,(1+tt)/4]; ap=[mp.mpf('0.5'),mp.mpf('-0.75'),mp.mpf('0.25')]
LA=[mp.log(tt)-mp.log(2),mp.log(1-tt)-2*mp.log(2),mp.log(1+tt)-2*mp.log(2)]
y=[(1-tt)*1+tt*2,(1-tt)*0+tt*2]; yp=[1,2]
E3_series = series_engine(tt,a,ap,LA,y,yp)
print("series engine E3(t=0.3):", mp.nstr(E3_series,30))

# ---- (2) verifier's symbolic t5 edge_E3 (K1c step_27 verbatim) ----
t,d = sp.symbols('t d', positive=True)
A = sp.symbols('A0:3'); X,Y,W = sp.symbols('X Y W'); l2=sp.log(2)
def edge_E3(V0,V1,terms):
    V0=sp.Matrix(V0); V1=sp.Matrix(V1)
    yy=(1-t)*V0+t*V1; ypp=V1-V0
    F=sum(aa*sp.exp(d*Aj) for (aa,Aj) in zip([q[0] for q in terms],A))
    Fp=(1+d)*sum(sp.diff(aa,t)*sp.exp(d*Aj) for (aa,Aj) in zip([q[0] for q in terms],A))
    Phi=Fp/((1+d)*F)
    G=F**(-2/(1+d))*(ypp.dot(ypp)-2*yy.dot(ypp)*Phi+yy.dot(yy)*Phi**2)
    ser=sp.series(sp.sqrt(G),d,0,4).removeO()
    E3=sp.expand(ser.coeff(d,3))
    sub={A[j]:terms[j][1] for j in range(len(terms))}
    return sp.expand(E3.xreplace(sub))
terms1=[(t/2,X-l2),(3*(1-t)/4,Y-2*l2),((1+t)/4,W-2*l2)]
E3_t5 = edge_E3((1,0),(2,2),terms1)
val_t5 = mp.mpf(sp.N(E3_t5.subs({t:tt,X:mp.log(tt),Y:mp.log(1-tt),W:mp.log(1+tt)},simultaneous=True),60))
print("t5 symbolic   E3(t=0.3):", mp.nstr(val_t5,30))

# ---- (3) FIXED pickle monomial evaluation, unweighted logs ----
E3pk = pickle.load(open('E3_fixed.pkl','rb'))
symtab = {s.name: s for s in sp.symbols('A0 A1 B L2 V02 kap lA0 lA1 lB q0')}
sub = {symtab['A0']:sp.Rational(1,2), symtab['A1']:sp.Rational(3,4), symtab['B']:sp.Rational(1,4),
       symtab['kap']:sp.Integer(1), symtab['L2']:sp.Integer(5), symtab['q0']:sp.Integer(1), symtab['V02']:sp.Integer(1),
       symtab['lA0']:sp.log(sp.Rational(1,2)), symtab['lA1']:sp.log(sp.Rational(1,4)), symtab['lB']:sp.log(sp.Rational(1,4))}
ttv=mp.mpf('0.3')
val_pk = mp.mpf(0)
for (n,a,b,c),co in E3pk['nl']+E3pk['loc']:
    cf = co.subs(sub, simultaneous=True)
    val_pk += mp.mpf(sp.N(cf,60)) * ttv**n * mp.log(ttv)**a * mp.log(1-ttv)**b * mp.log(1+ttv)**c
print("E3_fixed.pkl  E3(t=0.3):", mp.nstr(val_pk,30))

verif = mp.mpf('-0.64720423662744655873985211888')
print("\nverifier target          : -0.64720423662744655873985211888")
print("series - verifier :", mp.nstr(E3_series-verif,3))
print("t5     - verifier :", mp.nstr(val_t5-verif,3))
print("pkl    - verifier :", mp.nstr(val_pk-verif,3))

# ---- (4) same check for the OLD buggy pickle for contrast ----
E3old = pickle.load(open('E3_general.pkl','rb'))
val_old = mp.mpf(0)
for (n,a,b,c),co in E3old['nl']+E3old['loc']:
    cf = co.subs(sub, simultaneous=True)
    val_old += mp.mpf(sp.N(cf,60)) * ttv**n * mp.log(ttv)**a * mp.log(1-ttv)**b * mp.log(1+ttv)**c
print("OLD buggy pkl E3(t=0.3):", mp.nstr(val_old,30), " diff:", mp.nstr(val_old-verif,3))
'''
open('step_41.py','w').write(step_41)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_41.py'],
                   capture_output=True, text=True, timeout=900)
print(r.stdout); print(r.stderr[-2000:] if r.stderr else "")
with open('progress.log','a') as f:
    f.write("step_41: checked E3_fixed.pkl integrand at t=0.3 edge1 b=0 vs verifier -0.64720423662744655873985211888.\n"+r.stdout+"\n")
