evalsrc = '''"""Fast evaluator for Chen words G(w;1) on [0,1], letters mpc (0 allowed, != 1).
All singular points at distance >= ~0.414 from [0,1]; integrands analytic on [0,1].
Uniform panels, Gauss-Legendre per panel; rest-values F(u,v;x) tabulated at outer nodes
(cumulative panel sums + sub-panel quadrature); rest (u,0) via shuffle.
word_batch evaluates many words sharing rest tables.
"""
import mpmath as mp
from bisect import bisect_right
from collections import defaultdict

class GEval:
    def __init__(self, dps, M=64, gln=5, subn=5):
        mp.mp.dps=dps; self.dps=dps
        from mpmath.calculus.quadrature import GaussLegendre
        gl=GaussLegendre(mp.mp)
        prec=int(dps*3.34)+40
        self.gl=sorted(gl.calc_nodes(gln, prec), key=lambda p:p[0])
        self.glsub=sorted(gl.calc_nodes(subn, prec), key=lambda p:p[0])
        self.M=M
        self.edges=[mp.mpf(m)/M for m in range(M+1)]
        self.pts=[]   # all outer nodes (x,w), panel by panel
        for m in range(M):
            h2=mp.mpf(1)/(2*M); mid=(2*m+1)*h2
            for gx,gw in self.gl:
                self.pts.append((mid+h2*gx, gw*h2))
        self.restF={}   # (u,v) -> list of F values at pts
        self.restT={}   # (u,v) -> total
    def _h0(self, v, t):
        if t==0: return -1/v
        return mp.log1p(-t/v)/t
    def _f(self, u, v, t):
        if u==0: return self._h0(v,t)
        return mp.log(1-t/v)/(t-u)
    def _build_rest(self, u, v):
        if u==0 and v==0: raise ValueError("rest (0,0)")
        if u!=0 and v==0:
            arr0=self._build_rest(mp.mpc(0), u)
            arr=[mp.log(x)*mp.log(1-x/u)-arr0[j] for j,(x,w) in enumerate(self.pts)]
            self.restF[(u,v)]=arr
            self.restT[(u,v)]=-self.restT[(mp.mpc(0),u)]
            return arr
        # panel-edge cumulative
        cums={0:mp.mpc(0)}
        edges=self.edges; M=self.M
        glsub=self.glsub
        # compute F at all outer points: for each panel: cum at left edge via full GL, then
        # F(x) = cum_m + subGL(edge_m, x)
        arr=[None]*len(self.pts)
        cum=mp.mpc(0)
        for m in range(M):
            e0=edges[m]; e1=edges[m+1]
            panel=0
            s=mp.mpc(0)
            base=idx0=m*len(self.gl)
            for j in range(len(self.gl)):
                x,w=self.pts[idx0+j]
                h2=(x-e0)/2; mid=(x+e0)/2
                Fx=cum + h2*sum(w2*self._f(u,v,mid+h2*gx) for gx,w2 in glsub)
                arr[idx0+j]=Fx
                s+=w*self._f(u,v,x)
            cum+=s
        self.restF[(u,v)]=arr
        self.restT[(u,v)]=cum
        return arr
    def rest_total(self,u,v):
        if (u,v) not in self.restT: self._build_rest(u,v)
        return self.restT[(u,v)]
    def word_batch(self, words):
        """words: list of tuples (mpc letters, len 2 or 3). returns dict word->value."""
        need=defaultdict(list)
        for w in words: need[(w[1],w[2])].append(w)
        out={}
        for (u,v),ws in need.items():
            if (u,v) not in self.restF: self._build_rest(u,v)
            arr=self.restF[(u,v)]
            for w in ws:
                s=w[0]
                if s==0:
                    out[w]=sum(ww*arr[j]/x for j,(x,ww) in enumerate(self.pts))
                else:
                    out[w]=sum(ww*arr[j]/(x-s) for j,(x,ww) in enumerate(self.pts))
        return out
    def word(self,w):
        if len(w)==1:
            s=w[0]; assert s!=0
            return mp.log(1-1/s)
        if len(w)==2:
            return self.rest_total(w[0],w[1])
        return self.word_batch([w])[w]
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
import time
t0=time.time()
ge=GEval(80, M=48, gln=4, subn=4)
mx=0
for w in [(s0,sa,sa),(s0,sa,smia),(s0,smia,sa),(s0,smb,sa),(s0,sa,smb),(s0,sib,smib),(sa,smb,sia),(sa,s0,sib),(smb,sib,smib),(sia,smb,s0)]:
    wn=tuple(fembc(s) for s in w)
    mx=max(mx,abs(ge.word(wn)-Grec(w)))
print("max resid:", mp.nstr(mx,4), " time:", round(time.time()-t0,1),"s, gl nodes:", len(ge.gl), len(ge.glsub))
