import time, sys, json
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=80
import corner_series as cs
cc=[mp.mpf(s) for s in json.load(open('mk_pi2_k55.json'))]

# Custom evaluator: tanh-sinh nodes on [0,1], series once per node
from mpmath import TanhSinh
TS=TanhSinh(mp)

def mk_all(a, N, deg=6):
    c=mp.cos(a)
    # get nodes/weights for [0,1] via tanh-sinh (transformed)
    # Use TS.calc_nodes on a dummy to get level nodes; simpler: use mp's quad nodes manually
    # We'll integrate F_k(g)=[u^k]F over [0,1]; use per-interval TS
    cache={}
    def series_at(g):
        if g not in cache:
            cache[g]=cs.integrand_series(g,c,N)
        return cache[g]
    # nodes via TS for interval [0,1]
    nodes=TS.calc_nodes(deg, 10)  # abscissas in [-1,1]? mpmath returns (x,w) in [-1,1]
    return nodes, series_at, cache

# Simpler: use mp.quad but share cache across k by calling a single function returning list
def mk_shared(a, N, pts):
    c=mp.cos(a); cache={}
    def Fk(g,k):
        if g not in cache: cache[g]=cs.integrand_series(g,c,N)
        return cache[g][k]
    # Do all k in one quad by integrating vector? mp.quad is scalar. Use manual TS.
    # Manual tanh-sinh on each subinterval
    M=[mp.mpf(0)]*(N+1)
    def ts_integrate(f, a_, b_, n=400):
        # tanh-sinh: x=(a+b)/2+(b-a)/2*tanh(pi/2 sinh t)
        h=mp.mpf('0.02'); half=(b_-a_)/2; mid=(a_+b_)/2
        S=mp.mpf(0)
        for i in range(-n,n+1):
            t=i*h
            s=mp.sinh(t); th=mp.tanh(mp.pi/2*s)
            x=mid+half*th
            # dx/dt = half * (pi/2) cosh(t) / cosh(pi/2 sinh t)^2
            dx=half*(mp.pi/2)*mp.cosh(t)/(mp.cosh(mp.pi/2*s)**2)
            S+=f(x)*dx
        return S*h
    for k in range(N+1):
        M[k]=mp.fsum(ts_integrate(lambda g: Fk(g,k), pts[i],pts[i+1]) for i in range(len(pts)-1))
    return M,len(cache)

t=time.time()
M55,nn=mk_shared(mp.pi/2, 55, [mp.mpf(0),mp.mpf(1)/64,mp.mpf(1)/8,mp.mpf(1)/2,mp.mpf(1)])
print("time:",round(time.time()-t,1),"s, nodes used:",nn)
for k in [40,50,55]:
    print(f"k={k}: custom={mp.nstr(M55[k],14)} cc={mp.nstr(cc[k],14)} err={mp.nstr(abs(M55[k]-cc[k]),3)}")
