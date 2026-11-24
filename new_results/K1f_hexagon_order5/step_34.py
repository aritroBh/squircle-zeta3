code = r'''
import sympy as sp, pickle, time, sys
WD = '<project>/research/08_theory/K1f_hexagon_order5'
log = open(f'{WD}/step_24.out','w')
def pr(*a):
    s=' '.join(str(x) for x in a); log.write(s+'\n'); log.flush()

A0,A1,B,kap = sp.symbols('A0 A1 B kap')
L2,q0,V02 = sp.symbols('L2 q0 V02')
lA0,lA1,lB = sp.symbols('lA0 lA1 lB')
t,X,Y,W = sp.symbols('t X Y W')
N = 6

def mul(a,b):
    return [sp.expand(sum(a[i]*b[k-i] for i in range(k+1))) for k in range(N)]
def inv(a):
    r=[sp.Integer(1)]+[0]*(N-1)
    for k in range(1,N):
        r[k]=sp.expand(-sum(a[i]*r[k-i] for i in range(1,k+1)))
    return r
def log1(a):
    x=[sp.Integer(0)]+a[1:]
    res=[sp.Integer(0)]*N; pw=[sp.Integer(1)]+[0]*(N-1)
    for j in range(1,N):
        pw=mul(pw,x); res=[sp.expand(res[k]+sp.Rational((-1)**(j+1),j)*pw[k]) for k in range(N)]
    return res
def exp0(a):
    res=[sp.Integer(1)]+[0]*(N-1); pw=[sp.Integer(1)]+[0]*(N-1)
    for j in range(1,N):
        pw=mul(pw,a); res=[sp.expand(res[k]+pw[k]/sp.factorial(j)) for k in range(N)]
    return res

a0=A0*t; a1=A1*(1-t); a2=B*(1+kap*t)
LA0=X+lA0; LA1=Y+lA1; LA2=W+lB
t0=time.time()
T=[sp.Integer(1)]+[sp.expand(a0*LA0**k+a1*LA1**k+a2*LA2**k) for k in range(1,N)]
M=[sp.Integer(0)]+[sp.expand(A0*LA0**k-A1*LA1**k+B*kap*LA2**k) for k in range(1,N)]
F=[sp.expand(T[k]/sp.factorial(k)) for k in range(N)]
Num=[sp.expand(M[k]/sp.factorial(k)) for k in range(N)]
Phi=mul(Num,inv(F))
oneoverp=[sp.Integer((-1)**k) for k in range(N)]
Fpow=exp0([sp.expand(-c) for c in mul(oneoverp,log1(F))])
yyp=q0+L2*t; yy=V02+2*q0*t+L2*t**2
Phi2=mul(Phi,Phi)
Mmet=[sp.Integer(0)]*N
for k in range(N):
    val = -2*yyp*Phi[k] + yy*Phi2[k]
    if k==0: val += L2
    Mmet[k]=sp.expand(val)
Xmet=[sp.Integer(0)]+[sp.expand(Mmet[k]/L2) for k in range(1,N)]
sq=[sp.Integer(1)]+[0]*(N-1); pw=[sp.Integer(1)]+[0]*(N-1)
for j in range(1,N):
    pw=mul(pw,Xmet); sq=[sp.expand(sq[k]+sp.binomial(sp.Rational(1,2),j)*pw[k]) for k in range(N)]
Smet=[sp.expand(sp.sqrt(L2)*c) for c in sq]
E_a=mul(Fpow,Smet); E_a=[sp.expand(c) for c in E_a]
pr("engine built, %.1f s"%(time.time()-t0))
# sanity: E[0] == sqrt(L2)*1
pr("E[0] =", sp.simplify(E_a[0]))
# split into monomial dicts in (t,X,Y,W)
keyset=(t,X,Y,W)
polys={}
for k in range(6):
    D={}
    for term in sp.Add.make_args(E_a[k]):
        pd=term.as_powers_dict()
        key=(pd.get(t,0),pd.get(X,0),pd.get(Y,0),pd.get(W,0))
        rest=[s**e for s,e in pd.items() if s not in keyset]
        co=sp.Mul(*rest) if rest else sp.Integer(1)
        if key in D: D[key]=D[key]+co
        else: D[key]=co
    D={kk:sp.expand(v) for kk,v in D.items()}
    polys[k]=D
    pr(f"E[{k}]: {len(D)} monomials, {time.time()-t0:.0f}s")
pickle.dump(polys,open(f'{WD}/E5_polys.pkl','wb'))
pr("saved E5_polys.pkl")
# verify E[3] vs K1d
K1D = '<project>/research/08_theory/K1d_engine_fixed'
e3 = pickle.load(open(f'{K1D}/E3_fixed.pkl','rb'))
k1d={}
for mon,co in e3['nl']+e3['loc']:
    k1d[tuple(mon)]=sp.expand(co)
allk=set(polys[3])|set(k1d); bad=0
for k in sorted(allk):
    if sp.simplify(polys[3].get(k,sp.Integer(0))-k1d.get(k,sp.Integer(0)))!=0: bad+=1
pr("E[3] vs K1d: monomials mine", len(polys[3]), "k1d", len(k1d), "mismatches", bad)
with open(f'{WD}/progress.log','a') as f:
    f.write("step_24: rebuilt engine + E5_polys.pkl (dict form); E[3] vs K1d mismatches=%d; sizes=%s\n" % (bad,{k:len(v) for k,v in polys.items()}))
pr("DONE")
'''
wd = '<project>/research/08_theory/K1f_hexagon_order5'
open(f'{wd}/step_24.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
subprocess.Popen([PY, f'{wd}/step_24.py'], cwd=wd, stdout=open(f'{wd}/step_24.log','w'), stderr=subprocess.STDOUT)
print("launched step_24 in background")
