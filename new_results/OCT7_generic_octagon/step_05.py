
import sympy as sp, pickle, mpmath as mp
mp.mp.dps = 65
import mixed_engine_gen as E

mons = pickle.load(open('E3_edge1.pkl','rb'))['mons']
t = sp.symbols('t')
lA0, lA1, lB1, lB2 = sp.symbols('lA0 lA1 lB1 lB2')
# constant logs: lA0 = log(1/9) = -2 l3 ; lA1 = lB1 = log(1/4) = -2 l2 ; lB2 = log(1/2) = -l2
subs_logs = {lA0: -2*E.l3, lA1: -2*E.l2, lB1: -2*E.l2, lB2: -E.l2}

Earcs = sp.Integer(0)
for (exps, coef) in mons:
    a,b,c1,c2 = exps[:4]
    eA0,eA1,eB1,eB2 = exps[4:]
    # coefficient polynomial in t
    cpoly = sp.Poly(coef, t)
    tot = sp.Integer(0)
    for (n_,), cn in cpoly.terms():
        tot += cn * E.INT(n_, a,b,c1,c2)
    fac = lA0**eA0 * lA1**eA1 * lB1**eB1 * lB2**eB2
    Earcs += sp.expand(tot) * fac
Earcs = sp.expand(Earcs.xreplace(subs_logs))
print("E_arc assembled; number of terms:", len(sp.Add.make_args(Earcs)))

# numeric evaluation
def Li21(a, c, N=600):
    a, c = mp.mpf(a), mp.mpf(c)
    return mp.fsum([a**n * c**m / (mp.mpf(n)**2 * m) for n in range(2, N) for m in range(1, n)])
L2 = lambda v: mp.polylog(2, mp.mpf(v)); L3 = lambda v: mp.polylog(3, mp.mpf(v))
NUM = {
 E.l2: mp.log(2), E.l3: mp.log(3), E.l5: mp.log(5), E.z2: mp.pi**2/6, E.z3: mp.zeta(3),
 E.L2x1: L2(-mp.mpf(1)/3), E.L3x1: L3(-mp.mpf(1)/3), E.L2x2: L2(-mp.mpf(1)/9), E.L3x2: L3(-mp.mpf(1)/9),
 E.L2b1: L2(mp.mpf(1)/4), E.L3b1: L3(mp.mpf(1)/4), E.L2b2: L2(mp.mpf(1)/10), E.L3b2: L3(mp.mpf(1)/10),
 E.L2mh: L2(-mp.mpf(1)/2), E.L2m2t: L2(-mp.mpf(2)/3), E.L2m3h: L2(-mp.mpf(3)/2), E.L3m3h: L3(-mp.mpf(3)/2),
 E.L2m2: L2(-mp.mpf(2)), E.L3m2: L3(-mp.mpf(2)), E.L2_3f: L2(mp.mpf(3)/5), E.L3_3f: L3(mp.mpf(3)/5),
 E.L2_2t: L2(mp.mpf(2)/3), E.L3_2t: L3(mp.mpf(2)/3),
 E.Lam1: Li21(-mp.mpf(1)/3, mp.mpf(1)/3), E.Lam2: Li21(-mp.mpf(1)/9, mp.mpf(3)),
 E.Lam3: Li21(mp.mpf(1)/4, mp.mpf(2)/5), E.Lam4: Li21(mp.mpf(1)/10, mp.mpf(5)/2),
}
def ev(e): return mp.mpf(sp.N(e.subs(list(NUM.items())), mp.mp.dps))
Earc_num = ev(Earcs)
print("E_arc (exact assembly) =", mp.nstr(Earc_num, 40))

# independent quadrature of int_0^1 E3 dt
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
A0=mp.mpf(1)/9; A1=mp.mpf(1)/4; B1=mp.mpf(1)/4; B2=mp.mpf(1)/2
k1m=mp.mpf(1)/3; k2m=mp.mpf(1)/9
V0m=[mp.mpf(1)/4, mp.mpf(0)]; V1m=[mp.mpf(1)/3, -mp.mpf(1)/9]
def E3_num(tt):
    y=[V0m[0]+tt*(V1m[0]-V0m[0]), V0m[1]+tt*(V1m[1]-V0m[1])]
    yp=[V1m[0]-V0m[0], V1m[1]-V0m[1]]
    aa=[A0*tt, A1*(1-tt), B1*(1+k1m*tt), B2*(1+k2m*tt)]
    aap=[A0, -A1, B1*k1m, B2*k2m]
    LA=[mp.log(v) for v in aa]
    return series_engine(tt, aa, aap, LA, y, yp)
I = mp.quad(E3_num, [0,1])
print("independent quad int E3 dt =", mp.nstr(I, 40))
print("diff:", mp.nstr(abs(Earc_num-I),3))

# coefficient vector over named constants
from collections import OrderedDict
vec = OrderedDict()
for s in E.NAMED:
    vec[str(s)] = sp.simplify(Earcs.coeff(s))
print("\ncoefficients of named constants in E_arc:")
for k,v in vec.items():
    if v != 0: print(f"  {k}: {v}")
# the rest (products of logs etc.) - also print the residual part independent of NAMED polylogs
rest = sp.expand(Earcs - sum(Earcs.coeff(s)*s for s in E.NAMED))
print("\nresidual (log-monomial part):", rest)
pickle.dump({'Earcs': Earcs}, open('Earc_exact.pkl','wb'))
with open("progress.log","a") as f:
    f.write(f"STEP 05: E_arc assembled exactly; assembly vs independent quad diff {mp.nstr(abs(Earc_num-I),3)}\n")
