evalsrc = '''"""Fast evaluator for Chen words G(w;1), w over letters {0,a,-b,+-ia,+-ib} (as mpc),
on [0,1]. All integrands analytic on (0,1] with only integrable endpoint behavior at 0.
Scheme: outer tanh-sinh grid; depth-2 rests computed cumulatively with per-panel Gauss-Legendre;
rest (u,0) via shuffle G(u,0;t)=log t*log(1-t/u)-G(0,u;t).
"""
import mpmath as mp
from q1s_exact import feq

def ts_rule(h, T):
    h=mp.mpf(h); T=mp.mpf(T)
    K=int(mp.floor(T/h)); xs=[]; ws=[]
    for k in range(-K,K+1):
        t=k*h; u=(mp.pi/2)*mp.sinh(t); ch=mp.cosh(u)
        if u>0:
            e=mp.exp(-2*u); x=1/(1+e)
            if x==1: continue
        else:
            e=mp.exp(2*u); x=e/(1+e)
            if x==0: continue
        w=h*(mp.pi/2)*mp.cosh(t)/(2*ch*ch)
        xs.append(x); ws.append(w)
    return xs,ws

class GEval:
    def __init__(self, dps, h='0.008', T='4.8', gln=14):
        mp.mp.dps=dps; self.dps=dps
        self.xs,self.ws=ts_rule(h,T)
        self.N=len(self.xs)
        from mpmath.calculus.quadrature import GaussLegendre
        gl=GaussLegendre(mp.mp)
        self.glx,self.glw=gl.calc_nodes(gln, dps)
        self.cache1={}
        self.cache2={}
    def d1(self, u, t):
        # log(1 - t/u); u mpc != 0
        return mp.log(1-t/u)
    def _h0(self, v, t):
        # log(1-t/v)/t smooth at 0
        if t==0: return -1/v
        x=-t/v
        return mp.log1p(x)/t
    def rest_on_grid(self, u, v):
        """u,v mpc letters (u may be 0); returns (arr, total): G(u,v; x_k) and G(u,v;1)."""
        key=(u,v)
        if key in self.cache2: return self.cache2[key]
        xs=self.xs; N=self.N
        if u==0 and v==0: raise ValueError("rest (0,0)")
        if u!=0 and v==0:
            # shuffle: G(u,0;t) = log t * log(1-t/u) - G(0,u;t)
            arr0, tot0 = self.rest_on_grid(mp.mpc(0), u)
            arr=[mp.log(x)*self.d1(u,x)-arr0[k] for k,x in enumerate(xs)]
            total=mp.log(1)*0  # at t=1: log1=0 -> total = -G(0,u;1)
            total=-tot0
            self.cache2[key]=(arr,total)
            return self.cache2[key]
        if u==0:
            f=lambda t: self._h0(v,t)
        else:
            f=lambda t: self.d1(v,t)/(t-u)
        arr=[]; acc=mp.mpc(0); prev=mp.mpf(0)
        glx=self.glx; glw=self.glw
        for k in range(N):
            x1=xs[k]
            if x1==prev:
                arr.append(acc); continue
            h2=(x1-prev)/2; mid=(x1+prev)/2
            s=mp.mpc(0)
            for j in range(len(glx)):
                s+=glw[j]*f(mid+h2*glx[j])
            acc+=h2*s
            arr.append(acc); prev=x1
        self.cache2[key]=(arr,acc)
        return self.cache2[key]
    def word(self, w):
        """w tuple of mpc letters -> G(w;1)."""
        if len(w)==1:
            s=w[0]
            assert s!=0
            return mp.log(1-1/s)
        if len(w)==2:
            return self.rest_on_grid(w[0],w[1])[1]
        s=w[0]
        arr,_=self.rest_on_grid(w[1],w[2])
        tot=mp.mpc(0)
        if s==0:
            for k in range(self.N):
                tot+=self.ws[k]*arr[k]/self.xs[k]
        else:
            for k in range(self.N):
                tot+=self.ws[k]*arr[k]/(self.xs[k]-s)
        return tot
'''
open('q3_eval.py','w').write(evalsrc)
print("ok")
