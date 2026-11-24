import time, json
from fractions import Fraction as Fr
from math import comb, factorial
import sympy as sp

t0=time.time()
KK = 7   # compute mu_{k,n} for k<=KK, n<=30
u, y, a, b, n = sp.symbols('u y a b n')
A = -(n - 1)*a + n*b
l = b - a
# E_n = exp(uA) (1 + y(e^{u l} - 1))^{1-2n}; expand in u to order KK (symbolic n; flat_exact.py code path)
logE = u*A + (1 - 2*n)*sp.log(1 + y*(sp.exp(u*l) - 1))
logE_ser = sp.series(logE, u, 0, KK + 1).removeO()
E_ser = sp.series(sp.exp(logE_ser), u, 0, KK + 1).removeO()
print('E series done', round(time.time()-t0,1),'s', flush=True)
Ek = [sp.Poly(sp.expand(E_ser.coeff(u, k)), y, a, b) for k in range(KK + 1)]
print('Ek polys done', round(time.time()-t0,1),'s', flush=True)

# exact arithmetic in Q[z2..z7]: dict {(e2..e7): Fraction}; tuple length KK-1 (z2..z7)
ZL = KK-1
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
def psi(m, N):
    """psi^{(m)}(N) (m>=0), exact in Q[zetas]; psi(N) = H_{N-1} - gamma, gamma cancels in differences."""
    if m==0: return const(H(1,N-1))
    c = (-1)**(m+1)*factorial(m)
    return padd(pscale(zeta(m+1), c), const(-c*H(m+1,N-1)))

def cumulant(r, m, p, q):
    out = pscale(psi(r+m-1, p+q), -1)
    if m == 0: out = padd(out, psi(r+m-1, p))
    if r == 0: out = padd(out, psi(r+m-1, q))
    return out

_memo={}
def Dmom(i, l, p, q):
    key=(i,l,p,q)
    if key in _memo: return _memo[key]
    if i==0 and l==0: res=const(1)
    else:
        res={}
        if i>0:
            for r in range(1,i+1):
                for m in range(0,l+1):
                    t=pmul(cumulant(r,m,p,q), Dmom(i-r,l-m,p,q))
                    res=padd(res, pscale(t, comb(i-1,r-1)*comb(l,m)))
        else:
            for m in range(1,l+1):
                t=pmul(cumulant(0,m,p,q), Dmom(0,l-m,p,q))
                res=padd(res, pscale(t, comb(l-1,m-1)))
    _memo[key]=res
    return res

def poch_ratio(nv, j):
    r=Fr(1)
    for t in range(j): r*=Fr(nv+t,2*nv+t)
    return r

def mu(k, nv):
    tot={}
    for (j,i,l),co in Ek[k].terms():
        cr=sp.Rational(sp.expand(co.subs(n,nv))); cval=Fr(int(cr.p),int(cr.q))
        tot=padd(tot, pscale(Dmom(i,l,nv+j,nv), cval*poch_ratio(nv,j)))
    return pscale(tot, Fr(-1,nv*(2*nv-1)))

table={}
for k in range(KK+1):
    table[k]={}
    for nv in range(1,31):
        table[k][nv]={','.join(map(str,e)):str(v) for e,v in mu(k,nv).items()}
print('table done', round(time.time()-t0,1),'s', flush=True)
print('mu_{7,1} =', table[7][1])
print('mu_{7,2} =', table[7][2])
json.dump({str(k):{str(nv):table[k][nv] for nv in table[k]} for k in table}, open('mu7_indep_N30.json','w'))
print('saved mu7_indep_N30.json')
open('progress.log','a').write(f'step03a: independent symbolic-n mu_{{k,n}} (k<=7,n<=30) via sympy series path in {time.time()-t0:.0f}s; saved mu7_indep_N30.json\n')
