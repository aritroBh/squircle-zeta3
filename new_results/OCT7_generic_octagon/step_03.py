
import sympy as sp, pickle
from sympy import Rational, symbols, Matrix

# ---- edge 1 data (exact, from geometry.pkl / step_01) ----
V0 = Matrix([Rational(1,4), 0]); V1 = Matrix([Rational(1,3), Rational(-1,9)])
A0 = Rational(1,9); A1 = Rational(1,4)
B1 = Rational(1,4); k1 = Rational(1,3)
B2 = Rational(1,2); k2 = Rational(1,9)

t, d = symbols('t d')
L0, L1, L2, L3 = symbols('L0 L1 L2 L3')
a  = [A0*t, A1*(1-t), B1*(1+k1*t), B2*(1+k2*t)]
ap = [sp.diff(ai, t) for ai in a]
Ls_ = [L0, L1, L2, L3]
yy = (1-t)*V0 + t*V1
ypp = V1 - V0

n = 3
def expser(Li): return [1, Li, Li**2/2, Li**3/6]
F = [sum(a[i]*expser(Ls_[i])[k] for i in range(4)) for k in range(n+1)]
S = [sum(ap[i]*expser(Ls_[i])[k] for i in range(4)) for k in range(n+1)]
assert sp.expand(F[0]-1)==0, "F[0] != 1"
assert sp.expand(S[0])==0, "S[0] != 0"
def pmul(P,Q):
    R=[sp.Integer(0)]*(n+1)
    for i in range(n+1):
        for j in range(n+1-i): R[i+j]=sp.expand(R[i+j]+P[i]*Q[j])
    return R
def pdiv(P,Q):
    R=[sp.Integer(0)]*(n+1)
    for i in range(n+1):
        R[i]=sp.expand((P[i]-sum(R[j]*Q[i-j] for j in range(i)))/Q[0])
    return R
def plog(P):
    Pp=[P[i]*i for i in range(1,n+1)]+[sp.Integer(0)]
    Q=pdiv(Pp,P); L=[sp.Integer(0)]
    for i in range(n): L.append(sp.expand(Q[i]/(i+1)))
    return L
def pexp(P):
    E=[sp.Integer(1)]+[sp.Integer(0)]*n
    for k in range(1,n+1): E[k]=sp.expand(sum(P[j]*j*E[k-j] for j in range(1,k+1))/k)
    return E
Phi = pdiv(S, F)
lg = plog(F)
inv = [1,-1,1,-1]
m2lg = [sp.expand(-2*x) for x in lg]
fac = pexp(pmul(m2lg, inv))
yp2 = ypp.dot(ypp); yyp = yy.dot(ypp); yy2 = yy.dot(yy)
Phisq = pmul(Phi,Phi)
Q = [sp.Integer(0)]*(n+1)
Q[0] = yp2
for k in range(1,n+1): Q[k] = sp.expand(-2*yyp*Phi[k])
for k in range(n+1): Q[k] = sp.expand(Q[k] + yy2*Phisq[k])
G = pmul(fac, Q)
sq=[sp.Integer(0)]*(n+1); sq[0]=sp.sqrt(G[0])
for k in range(1,n+1): sq[k]=sp.expand((G[k]-sum(sq[j]*sq[k-j] for j in range(1,k)))/(2*sq[0]))
E3 = sp.expand(sq[3])
print("E3 built; |y'| =", sp.sqrt(yp2), "; terms:", len(sp.Add.make_args(E3)))

# ---- substitute logs: L0=X+lA0, L1=Y+lA1, L2=W1+lB1, L3=W2+lB2 ----
X, Y, W1, W2 = symbols('X Y W1 W2')
lA0, lA1, lB1, lB2 = symbols('lA0 lA1 lB1 lB2')
E3s = sp.expand(E3.xreplace({L0: X+lA0, L1: Y+lA1, L2: W1+lB1, L3: W2+lB2}))
P = sp.Poly(E3s, X, Y, W1, W2, lA0, lA1, lB1, lB2)
mons = P.terms()
print("log-monomials (incl constant logs):", len(mons))
# classify
from collections import Counter
cls = Counter()
for (exps, coef) in mons:
    a_,b_,c1,c2 = exps[:4]
    kind = 'local' if (c1==0 and c2==0) else ('s1' if c2==0 else ('s2' if c1==0 else 'MIXED'))
    cls[kind]+=1
print(cls)

# ---- numeric verification at t=0.3, 0.7 against direct numeric engine ----
import mpmath as mp
mp.mp.dps = 50
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
    Phi=pdiv(Fp,F); lg=plog(F); inv=[mp.mpf(1),mp.mpf(-1),mp.mpf(1),mp.mpf(-1)][:n+1]
    fac=pexp(pmul([mp.mpf(-2)*x for x in lg],inv))
    yp2=yp[0]**2+yp[1]**2; yyp=y[0]*yp[0]+y[1]*yp[1]; yy2=y[0]**2+y[1]**2
    Q=[mp.mpf(0)]*(n+1); Q[0]=yp2
    for k in range(1,n+1): Q[k]=-2*yyp*Phi[k]
    Phisq=pmul(Phi,Phi)
    for k in range(n+1): Q[k]+=yy2*Phisq[k]
    G=pmul(fac,Q)
    sq=[mp.mpf(0)]*(n+1); sq[0]=mp.sqrt(G[0])
    for k in range(1,n+1): sq[k]=(G[k]-sum(sq[j]*sq[k-j] for j in range(1,k)))/(2*sq[0])
    return sq[3]
def E3_num(tt):
    y = [V0[0]+tt*(V1[0]-V0[0]), V0[1]+tt*(V1[1]-V0[1])]
    y = [mp.mpf(str(v)) for v in y]
    yp = [mp.mpf(str(V1[0]-V0[0])), mp.mpf(str(V1[1]-V0[1]))]
    aa = [A0*tt, A1*(1-tt), B1*(1+k1*tt), B2*(1+k2*tt)]
    aa = [mp.mpf(str(v)) for v in aa]
    aap = [A0, -A1, B1*k1, B2*k2]; aap=[mp.mpf(str(v)) for v in aap]
    LA = [mp.log(v) for v in aa]
    return series_engine(tt, aa, aap, LA, y, yp)
vals = {lA0: mp.log(mp.mpf(1)/9), lA1: mp.log(mp.mpf(1)/4), lB1: mp.log(mp.mpf(1)/4), lB2: mp.log(mp.mpf(1)/2)}
worst = mp.mpf(0)
for ttq in [mp.mpf('0.2'), mp.mpf('0.3'), mp.mpf('0.55'), mp.mpf('0.8')]:
    tt = mp.mpf(str(ttq))
    sub = {t: tt, X: mp.log(tt), Y: mp.log(1-tt),
           W1: mp.log(1+mp.mpf(1)/3*tt), W2: mp.log(1+mp.mpf(1)/9*tt)}
    sub.update(vals)
    vsym = mp.mpf(sp.N(E3s.xreplace(sub), 50))
    vnum = E3_num(tt)
    dff = abs(vsym-vnum); worst = max(worst, dff)
    print(f"t={tt}: sym={mp.nstr(vsym,25)} num={mp.nstr(vnum,25)} diff={mp.nstr(dff,3)}")
print("worst E3 point-check diff:", mp.nstr(worst,3))
pickle.dump({'mons': mons, 'E3s': E3s}, open('E3_edge1.pkl','wb'))
with open("progress.log","a") as f:
    f.write(f"STEP 03: exact E3 for edge 1 built ({len(mons)} monomials, |y'|=5/36); point checks worst {mp.nstr(worst,3)}\n")
print("saved E3_edge1.pkl")
