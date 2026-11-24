
import mpmath as mp, json, time
WORK='research/08_theory/K1e_hexagon_closed'
mp.mp.dps = 100
t0=time.time()
NJ=4
def madd(a,b): return [x+y for x,y in zip(a,b)]
def mscal(c,a): return [c*x for x in a]
def mmul(a,b): return [mp.fsum([a[i]*b[k-i] for i in range(k+1)]) for k in range(NJ)]
def minv(a):
    r=[1/a[0]]+[mp.mpf(0)]*3
    for k in range(1,NJ): r[k]=-mp.fsum([a[i]*r[k-i] for i in range(1,k+1)])/a[0]
    return r
def mlog(a):
    x=[mp.mpf(0)]+[a[k]/a[0] for k in range(1,NJ)]
    res=[mp.mpf(0)]*NJ; pw=[mp.mpf(1)]+[mp.mpf(0)]*3
    for j in range(1,NJ):
        pw=mmul(pw,x); res=[res[k]+mp.mpf((-1)**(j+1))/j*pw[k] for k in range(NJ)]
    res[0]=mp.log(a[0]); return res
def mexp(a):
    x=[mp.mpf(0)]+list(a[1:])
    res=[mp.mpf(1)]+[mp.mpf(0)]*3; pw=list(res)
    for j in range(1,NJ):
        pw=mmul(pw,x); res=[res[k]+pw[k]/mp.factorial(j) for k in range(NJ)]
    e0=mp.exp(a[0]); return [e0*r for r in res]
def msqrt(a): return mexp(mscal(mp.mpf('0.5'), mlog(a)))

# polar jet integrand from the DEFINITION: rho = F^(-1/p), p=1+d; sqrt(rho^2+rho'^2)
vs = [(mp.mpf(0), mp.mpf(1)/18), (-mp.mpf(3)/18, mp.mpf(3)/18), (-mp.mpf(1)/18, mp.mpf(0))]
cs = [mp.mpf(3), mp.mpf(5), mp.mpf(3)]
def polar_jet(th):
    u=(mp.cos(th), mp.sin(th)); du=(-mp.sin(th), mp.cos(th))
    Fj=[mp.mpf(0)]*NJ; dFj=[mp.mpf(0)]*NJ
    ip=[mp.mpf(1), -mp.mpf(1), mp.mpf(1), -mp.mpf(1)]  # 1/(1+d)
    for v,c in zip(vs,cs):
        g=v[0]*u[0]+v[1]*u[1]; gp=v[0]*du[0]+v[1]*du[1]
        s=abs(g); sg=1 if g>0 else -1
        ls=mp.log(s)
        term=mexp([mp.log(c*s), ls, mp.mpf(0), mp.mpf(0)])   # c*s*exp(d*log s)
        Fj=madd(Fj, term)
        t2=mexp([mp.log(c), ls, mp.mpf(0), mp.mpf(0)])
        t2=mmul(t2, [mp.mpf(1), mp.mpf(1), mp.mpf(0), mp.mpf(0)])  # *(1+d)
        dFj=madd(dFj, mscal(sg*gp, t2))
    lF=mlog(Fj)
    rho=mexp(mscal(-1, mmul(ip, lF)))
    drho=mscal(-1, mmul(mmul(rho, ip), mmul(dFj, minv(Fj))))
    return msqrt(madd(mmul(rho,rho), mmul(drho,drho)))

import sympy as sp
geo = json.load(open(f"{WORK}/geometry_fresh.json"))
direct={}
for e in geo:
    V0x, V0y = [mp.mpf(sp.Integer(c)) for c in e['V0']]
    V1x, V1y = [mp.mpf(sp.Integer(c)) for c in e['V1']]
    Dx, Dy = V1x-V0x, V1y-V0y
    cr = V0x*Dy - V0y*Dx   # cross(y,D) = cross(V0,D) (t-independent!)
    def g(t, cr=cr, V0x=V0x, V0y=V0y, Dx=Dx, Dy=Dy):
        yx = V0x + t*Dx; yy = V0y + t*Dy
        th = mp.atan2(yy, yx)
        dth = cr/(yx*yx+yy*yy)
        return polar_jet(th)[3]*dth
    I = mp.quad(g, [0,1], method='tanh-sinh')
    P3 = 6*I
    direct[e['edge']] = P3
    print(f"edge {e['edge']} (r={e['r']}): direct P3_arc = {mp.nstr(P3,45)}  ({time.time()-t0:.0f}s)")
tot = 2*sum(direct.values())
print("direct total P'''(1) =", mp.nstr(tot,45))
json.dump({str(k): mp.nstr(v,95) for k,v in direct.items()}, open(f"{WORK}/fresh_direct.json","w"))
