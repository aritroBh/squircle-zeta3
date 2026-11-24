import time, sys, json
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
import corner_series as cs

def mk_shared(a, N, pts, h, n):
    c=mp.cos(a); cache={}
    def Fk(g,k):
        if g not in cache: cache[g]=cs.integrand_series(g,c,N)
        return cache[g][k]
    tiny=mp.mpf(10)**(-mp.mp.dps+10)
    def ts_integrate(f, a_, b_):
        half=(b_-a_)/2; mid=(a_+b_)/2; S=mp.mpf(0)
        for i in range(-n,n+1):
            t=i*h; s=mp.sinh(t); th=mp.tanh(mp.pi/2*s)
            x=mid+half*th
            if x<=a_+tiny or x>=b_-tiny: continue
            dx=half*(mp.pi/2)*mp.cosh(t)/(mp.cosh(mp.pi/2*s)**2)
            S+=f(x)*dx
        return S*h
    M=[mp.mpf(0)]*(N+1)
    for k in range(N+1):
        M[k]=mp.fsum(ts_integrate(lambda g: Fk(g,k), pts[i],pts[i+1]) for i in range(len(pts)-1))
    return M,len(cache)

pts=[mp.mpf(0),mp.mpf(1)/64,mp.mpf(1)/8,mp.mpf(1)/2,mp.mpf(1)]
for dps in [100,140]:
    mp.mp.dps=dps
    t=time.time()
    M,nn=mk_shared(mp.pi/2,100,pts,mp.mpf('0.015'),350)
    print(f"dps={dps}: time={round(time.time()-t,1)}s nodes={nn}")
    for k in [60,80,100]:
        print(f"   M_{k}={mp.nstr(M[k],20)}")
    json.dump([mp.nstr(x,60) for x in M], open(f'mk_pi2_k100_dps{dps}.json','w'))
