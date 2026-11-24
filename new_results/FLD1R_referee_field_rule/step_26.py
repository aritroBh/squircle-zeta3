import time, json
from fractions import Fraction as Fr
from math import comb, factorial
import sympy as sp

t0=time.time()
KK=7; NMAX=30
u,y,a,b = sp.symbols('u y a b')

# ---- per-n exact E-series (independent code path: sympy series with numeric n)
def Ek_at(k, nv):
    A = -(nv - 1)*a + nv*b
    l = b - a
    logE = u*A + (1 - 2*nv)*sp.log(1 + y*(sp.exp(u*l) - 1))
    E_ser = sp.series(sp.exp(logE), u, 0, k+1).removeO()
    return sp.Poly(sp.expand(E_ser.coeff(u, k)), y, a, b)

ZL=KK-1
def padd(P,Q,s=1):
    R=dict(P)
    for k,v in Q.items():
        R[k]=R.get(k,0)+s*v
        if R[k]==0: del R[k]
    return R
def pmul(P,Q):
    R={}
    for k1,v1 in P.items():
        for k2,v2 in Q.items():
            kk=tuple(x+y_ for x,y_ in zip(k1,k2)); R[kk]=R.get(kk,0)+v1*v2
            if R[kk]==0: del R[kk]
    return R
def pscale(P,c): return {k:v*c for k,v in P.items()} if c else {}
def const(x): return {tuple([0]*ZL): Fr(x)} if x else {}
def zeta(idx):
    e=[0]*ZL; e[idx-2]=1; return {tuple(e): Fr(1)}
_H={}
def H(r,m):
    key=(r,m)
    if key not in _H: _H[key]=sum((Fr(1,i**r) for i in range(1,m+1)),Fr(0))
    return _H[key]
def psi(m,N):
    if m==0: return const(H(1,N-1))
    c=(-1)**(m+1)*factorial(m)
    return padd(pscale(zeta(m+1),c), const(-c*H(m+1,N-1)))
def cumulant(r,m,p,q):
    out=pscale(psi(r+m-1,p+q),-1)
    if m==0: out=padd(out,psi(r+m-1,p))
    if r==0: out=padd(out,psi(r+m-1,q))
    return out
_memo={}
def Dmom(i,l,p,q):
    key=(i,l,p,q)
    if key in _memo: return _memo[key]
    if i==0 and l==0: res=const(1)
    else:
        res={}
        if i>0:
            for r in range(1,i+1):
                for m in range(0,l+1):
                    res=padd(res, pscale(pmul(cumulant(r,m,p,q),Dmom(i-r,l-m,p,q)), comb(i-1,r-1)*comb(l,m)))
        else:
            for m in range(1,l+1):
                res=padd(res, pscale(pmul(cumulant(0,m,p,q),Dmom(0,l-m,p,q)), comb(l-1,m-1)))
    _memo[key]=res
    return res

def mu(k,nv):
    tot={}
    Ek=Ek_at(k,nv)
    for (j,i,l),co in Ek.terms():
        cval=Fr(int(co.p),int(co.q))
        pr=Fr(1)
        for t in range(j): pr*=Fr(nv+t,2*nv+t)
        tot=padd(tot, pscale(Dmom(i,l,nv+j,nv), cval*pr))
    return pscale(tot, Fr(-1,nv*(2*nv-1)))

table={}
for k in range(KK+1):
    table[k]={}
    for nv in range(1,NMAX+1):
        table[k][nv]={','.join(map(str,e)):str(v) for e,v in mu(k,nv).items()}
    print(f'k={k} done ({time.time()-t0:.0f}s)', flush=True)
print('mu_{7,1} =', table[7][1])
print('mu_{7,3} =', table[7][3])
json.dump({str(k):{str(nv):table[k][nv] for nv in table[k]} for k in table}, open('mu7_indep_N30.json','w'))
print('saved mu7_indep_N30.json in', round(time.time()-t0,1),'s')
