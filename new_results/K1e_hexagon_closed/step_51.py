
import mpmath as mp, json, time
WORK='research/08_theory/K1e_hexagon_closed'
mp.mp.dps = 100
t0=time.time()
NJ=4
def madd(a,b): return [x+y for x,y in zip(a,b)]
def msub(a,b): return [x-y for x,y in zip(a,b)]
def mscal(c,a): return [c*x for x in a]
def mmul(a,b): return [mp.fsum([a[i]*b[k-i] for i in range(k+1)]) for k in range(NJ)]
def minv(a):
    r=[1/a[0]]+[mp.mpf(0)]*3
    for k in range(1,NJ):
        r[k]=-mp.fsum([a[i]*r[k-i] for i in range(1,k+1)])/a[0]
    return r
def mlog(a):
    x=[mp.mpf(0)]+[a[k]/a[0] for k in range(1,NJ)]
    res=[mp.mpf(0)]*NJ; pw=[mp.mpf(1)]+[mp.mpf(0)]*3
    for j in range(1,NJ):
        pw=mmul(pw,x)
        res=[res[k]+mp.mpf((-1)**(j+1))/j*pw[k] for k in range(NJ)]
    res[0]=mp.log(a[0]); return res
def mexp(a):
    x=[mp.mpf(0)]+list(a[1:])
    res=[mp.mpf(1)]+[mp.mpf(0)]*3; pw=list(res)
    for j in range(1,NJ):
        pw=mmul(pw,x)
        res=[res[k]+pw[k]/mp.factorial(j) for k in range(NJ)]
    e0=mp.exp(a[0]); return [e0*r for r in res]
def msqrt(a): return mexp(mscal(mp.mpf('0.5'), mlog(a)))

import sympy as sp
geo = json.load(open(f"{WORK}/geometry_fresh.json"))
def edge_c3_integrand(e):
    A0=mp.mpf(sp.Rational(e['A0'])); A1=mp.mpf(sp.Rational(e['A1'])); B=mp.mpf(sp.Rational(e['B']))
    kap=mp.mpf(sp.Rational(e['kap'])); L2=mp.mpf(sp.Integer(e['L2'])); q0=mp.mpf(sp.Integer(e['q0'])); V02=mp.mpf(sp.Integer(e['V02']))
    w=[mp.mpf(sp.Rational(e['w0'])), mp.mpf(sp.Rational(e['w1'])), mp.mpf(sp.Rational(e['wo']))]
    invp=[mp.mpf(1)/2, -mp.mpf(1)/4, mp.mpf(1)/8, -mp.mpf(1)/16]
    m2=mscal(-2, invp)
    def g(t):
        av=[A0*t, A1*(1-t), B*(1+kap*t)]
        uv=[av[j]/w[j] for j in range(3)]
        LA=[mp.log(u) for u in uv]
        F=[mp.fsum([av[j]*uv[j]*LA[j]**k/mp.factorial(k) for j in range(3)]) for k in range(NJ)]
        Fp=[(k+1)*F[k+1] if k<3 else mp.mpf(0) for k in range(NJ)]
        Phi=mmul(Fp, mmul(invp, minv(F)))
        yD=q0+t*L2; yy=V02+2*q0*t+L2*t*t
        S=msub(madd([L2,0,0,0], mscal(-2*yD, Phi)), mscal(-yy, mmul(Phi,Phi)))
        Fpow=mexp(mmul(m2, mlog(F)))
        G=mmul(Fpow, S)
        return msqrt(G)[3]
    return g
results={}
for e in geo:
    g=edge_c3_integrand(e)
    I=mp.quad(g, [0,1], method='tanh-sinh')
    P3=6*I
    results[e['edge']]=P3
    print(f"edge {e['edge']} (r={e['r']}): direct P3_arc = {mp.nstr(P3,50)}   ({time.time()-t0:.0f}s)")
tot=2*sum(results.values())
print("direct P'''(1) total =", mp.nstr(tot,50))
json.dump({str(k): mp.nstr(v,95) for k,v in results.items()}, open(f"{WORK}/fresh_direct.json","w"))
print(f"({time.time()-t0:.0f}s)")
