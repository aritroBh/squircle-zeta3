evalsrc = '''"""Fast evaluator for Chen words G(w;1) on [0,1], letters mpc (0 allowed).
Integrands analytic on [0,1] (singular points at distance >= ~0.414 off-path).
Scheme: uniform panels with high-order Gauss-Legendre.
- rest(u,v) returns (F, total) with F(x) = G(u,v;x) computed as
  cum[panel(x)] + GL-quadrature of f on [edge, x]  (panel edges exact-ish).
- rest (u,0) via shuffle G(u,0;t) = log t*log(1-t/u) - G(0,u;t).
"""
import mpmath as mp
from bisect import bisect_right

class GEval:
    def __init__(self, dps, M=64, gln=4, subn=3):
        mp.mp.dps=dps; self.dps=dps
        from mpmath.calculus.quadrature import GaussLegendre
        gl=GaussLegendre(mp.mp)
        prec=int(dps*3.34)+40
        self.gl=sorted(gl.calc_nodes(gln, prec), key=lambda p:p[0])   # 3*2^gln nodes
        self.glsub=sorted(gl.calc_nodes(subn, prec), key=lambda p:p[0])
        self.M=M
        self.edges=[mp.mpf(m)/M for m in range(M+1)]
        # global GL points for panel sums
        self.pan=[]   # per panel: list of (x,w)
        for m in range(M):
            h2=(self.edges[m+1]-self.edges[m])/2; mid=(self.edges[m+1]+self.edges[m])/2
            self.pan.append([(mid+h2*gx, gw*h2) for gx,gw in self.gl])
        self.cache2={}
    def _h0(self, v, t):
        if t==0: return -1/v
        return mp.log1p(-t/v)/t
    def _quad_gl(self, f, x0, x1, nodes):
        h2=(x1-x0)/2; mid=(x1+x0)/2
        return h2*sum(w*f(mid+h2*gx) for gx,w in nodes)
    def rest(self, u, v):
        """returns (F, total): F(x)=G(u,v;x), total=G(u,v;1)."""
        key=(u,v)
        if key in self.cache2: return self.cache2[key]
        if u==0 and v==0: raise ValueError("rest (0,0)")
        if u!=0 and v==0:
            F0, tot0 = self.rest(mp.mpc(0), u)
            F=lambda x: mp.log(x)*mp.log(1-x/u)-F0(x)
            res=(F, -tot0)
            self.cache2[key]=res
            return res
        if u==0:
            f=lambda t: self._h0(v,t)
        else:
            f=lambda t: mp.log(1-t/v)/(t-u)
        # panel-edge cumulative sums
        cum=[mp.mpc(0)]
        for m in range(self.M):
            s=sum(w*f(x) for x,w in self.pan[m])
            cum.append(cum[-1]+s)
        total=cum[-1]
        edges=self.edges; M=self.M; glsub=self.glsub; _q=self._quad_gl
        def F(x):
            m=bisect_right(edges, x)-1
            if m>=M: m=M-1
            e=edges[m]
            return cum[m] + _q(f, e, x, glsub)
        res=(F,total)
        self.cache2[key]=res
        return res
    def word(self, w):
        if len(w)==1:
            s=w[0]; assert s!=0
            return mp.log(1-1/s)
        if len(w)==2:
            return self.rest(w[0],w[1])[1]
        s=w[0]
        F,_=self.rest(w[1],w[2])
        tot=mp.mpc(0)
        for m in range(self.M):
            if s==0:
                tot+=sum(ww*F(x)/x for x,ww in self.pan[m])
            else:
                tot+=sum(ww*F(x)/(x-s) for x,ww in self.pan[m])
        return tot
'''
open('q3_eval.py','w').write(evalsrc)
import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp
mp.mp.dps=80
def fembc(x): return femb(x,1)
def Grec(w, z=1):
    if len(w)==1:
        s=fembc(w[0])
        if feq(w[0],ZERO): return mp.log(z)
        return mp.log(1-z/s)
    s=fembc(w[0])
    return mp.quad(lambda t: Grec(w[1:],t)/(t-s), [0,z])
ge=GEval(80, M=48, gln=4, subn=3)
for w in [(s0,sa,sa),(s0,sa,smia),(s0,smia,sa),(s0,smb,sa),(s0,sa,smb),(s0,sib,smib),(sa,smb,sia),(sa,s0,sib)]:
    wn=tuple(fembc(s) for s in w)
    v1=ge.word(wn); v0=Grec(w)
    print(mp.nstr(abs(v1-v0),4), " ", end="")
print()
