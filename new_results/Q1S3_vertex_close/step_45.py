evalsrc = '''"""Fast evaluator for Chen words G(w;1) on [0,1], letters mpc (0 allowed).
All integrands analytic on [0,1] (singular points at distance >= ~0.414 off-path).
Scheme: uniform panels with high-order Gauss-Legendre; depth-2 rests cached at all GL nodes;
rest (u,0) via shuffle G(u,0;t) = log t*log(1-t/u) - G(0,u;t).
"""
import mpmath as mp

class GEval:
    def __init__(self, dps, M=64, gln=4):
        mp.mp.dps=dps; self.dps=dps
        from mpmath.calculus.quadrature import GaussLegendre
        gl=GaussLegendre(mp.mp)
        self.gl=gl.calc_nodes(gln, int(dps*3.34)+40)   # 3*2^gln nodes on [-1,1]
        self.M=M
        # panel edges and GL points
        n=len(self.gl)
        self.pts=[]   # (x, w) global GL points
        for m in range(M):
            h2=mp.mpf(1)/(2*M); mid=(2*m+1)*h2
            for gx,gw in self.gl:
                self.pts.append((mid+h2*gx, gw*h2))
        self.N=len(self.pts)
        self.cache2={}
    def _h0(self, v, t):
        if t==0: return -1/v
        return mp.log1p(-t/v)/t
    def rest_all(self, u, v):
        """returns (arr, total): G(u,v; x_j) at all self.pts, and G(u,v;1)."""
        key=(u,v)
        if key in self.cache2: return self.cache2[key]
        if u==0 and v==0: raise ValueError("rest (0,0)")
        if u!=0 and v==0:
            arr0, tot0 = self.rest_all(mp.mpc(0), u)
            arr=[mp.log(x)*mp.log(1-x/u)-arr0[j] for j,(x,w) in enumerate(self.pts)]
            res=(arr,-tot0)
            self.cache2[key]=res
            return res
        if u==0:
            f=lambda t: self._h0(v,t)
        else:
            f=lambda t: mp.log(1-t/v)/(t-u)
        n=len(self.gl); M=self.M
        arr=[None]*(M*n); cum=mp.mpc(0); idx=0
        for m in range(M):
            base=cum; vals=[]
            for j in range(n):
                x,w=self.pts[idx]
                fv=f(x)
                vals.append((w,fv))
                arr[idx]=base+sum(vv[0]*vv[1] for vv in vals)-vals[-1][0]*vals[-1][1]  # cum + partial
                idx+=1
            # recompute node values exactly: arr[idx] should be base + sum_{j'<=j} w f
            # do it cleanly:
            base2=base; s=mp.mpc(0)
            for j in range(n):
                s+=vals[j][0]*vals[j][1]
                arr[m*n+j]=base2+s
            cum=base2+s
        self.cache2[key]=(arr,cum)
        return self.cache2[key]
    def word(self, w):
        if len(w)==1:
            s=w[0]; assert s!=0
            return mp.log(1-1/s)
        if len(w)==2:
            return self.rest_all(w[0],w[1])[1]
        s=w[0]
        arr,_=self.rest_all(w[1],w[2])
        tot=mp.mpc(0)
        if s==0:
            for j,(x,ww) in enumerate(self.pts):
                tot+=ww*arr[j]/x
        else:
            for j,(x,ww) in enumerate(self.pts):
                tot+=ww*arr[j]/(x-s)
        return tot
'''
open('q3_eval.py','w').write(evalsrc)
import importlib, q3_eval
importlib.reload(q3_eval)
import mpmath as mp
mp.mp.dps=130
r2=mp.sqrt(2)
cases = []
pairs = [(1j*(r2-1),'ib'),(-1j*(r2-1),'-ib'),(1+r2,'a'),(-(r2-1),'-b'),(1j*(1+r2),'ia'),(-1j*(1+r2),'-ia')]
for v,name in pairs:
    cases.append(((mp.mpc(0),mp.mpc(v)), -mp.polylog(2,1/v)))
    cases.append(((mp.mpc(v),mp.mpc(v)), mp.log(1-1/v)**2/2))
for M,gln in [(32,4),(64,4),(64,5)]:
    ge=q3_eval.GEval(130,M=M,gln=gln)
    mx=0
    for w,val in cases:
        mx=max(mx,abs(ge.word(w)-val))
    print("M=%d gln=%d (n=%d pts): max resid %s"%(M,gln,ge.N, mp.nstr(mx,4)))
