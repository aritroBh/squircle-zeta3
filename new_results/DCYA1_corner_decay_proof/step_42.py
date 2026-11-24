import time, sys, json
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=80
import corner_series as cs
cc=[mp.mpf(s) for s in json.load(open('mk_pi2_k55.json'))]

def mk_shared(a, N, pts, h=mp.mpf('0.015'), n=350):
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

t=time.time()
M55,nn=mk_shared(mp.pi/2, 55, [mp.mpf(0),mp.mpf(1)/64,mp.mpf(1)/8,mp.mpf(1)/2,mp.mpf(1)])
print("time:",round(time.time()-t,1),"s, nodes used:",nn)
for k in [30,40,50,55]:
    print(f"k={k}: custom={mp.nstr(M55[k],16)} cc={mp.nstr(cc[k],16)} err={mp.nstr(abs(M55[k]-cc[k]),3)}")
