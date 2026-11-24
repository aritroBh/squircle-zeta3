import sympy as sp, mpmath as mp
import pickle
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
with open(f'{WD}/golden_edges.pkl','rb') as f: GD = pickle.load(f)
k_sym=GD['k']; V=GD['V']; w=GD['w']
mp.mp.dps=80
kk=[mp.mpf(sp.N(kj,90)) for kj in k_sym]
Vn=[(mp.mpf(sp.N(v[0],90)),mp.mpf(sp.N(v[1],90))) for v in V]
def P(p):
    cc=[w[j]*kk[j]**p for j in range(3)]
    def g(th):
        s,co=mp.sin(th),mp.cos(th); F=mp.mpf(0); Fp=mp.mpf(0)
        for j in range(3):
            sj=Vn[j][0]*s-Vn[j][1]*co; dsj=Vn[j][0]*co+Vn[j][1]*s
            F+=cc[j]*abs(sj)**p; Fp+=cc[j]*p*abs(sj)**(p-1)*mp.sign(sj)*dsj
        rho=F**(-mp.mpf(1)/p); drho=-(mp.mpf(1)/p)*F**(-mp.mpf(1)/p-1)*Fp
        return mp.sqrt(rho**2+drho**2)
    br=[mp.atan2(v[1],v[0]) for v in Vn]+[mp.pi]
    return 2*sum(mp.quad(g,[br[i],br[i+1]]) for i in range(3))

# 5th derivative weights for antisymmetric stencil points h,2h,3h,4h
mm=sp.symbols('m1:5')
eqs=[sp.Eq(sum(mm[i]*(i+1)**p for i in range(4)), val) for p,val in [(1,0),(3,0),(5,60),(7,0)]]
wsol=sp.solve(eqs,mm,dict=True)[0]
cw=[sp.nsimplify(wsol[mm[i]]) for i in range(4)]
print("5th-deriv antisymmetric weights (k=1..4):", cw)
def d5(h):
    return sum(cw[i]*(P(1+(i+1)*h)-P(1-(i+1)*h)) for i in range(4))/h**5
# Richardson
for h0 in [mp.mpf(10)**-4, mp.mpf(10)**-5]:
    a=d5(h0); b=d5(h0/2)
    # base accuracy O(h^4): richardson (16b-a)/15
    rich=(16*b-a)/15
    print(f"h0={h0}: d5(h0)={mp.nstr(a,18)}, d5(h0/2)={mp.nstr(b,18)}, Richardson={mp.nstr(rich,20)}")
