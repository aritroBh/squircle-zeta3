import mpmath as mp
mp.mp.dps=80
import sympy as sp, pickle
from poly_coeffs import regular_zeros
import sym_G
z3=mp.zeta(3); pi=mp.pi
# 1) plain single-corner log^3 sine integral
I=mp.quad(lambda th: mp.log(mp.sin(th))**3,[0,mp.pi/2],method='tanh-sinh')
print('int_0^{pi/2} log^3(sin th) dth =',mp.nstr(I,30))
basis=[z3,pi**3,pi**2*mp.log(2),pi*mp.log(2)**2,mp.log(2)**3,pi*mp.catalan]
bn=['z3','pi^3','pi^2log2','pi*log2^2','log2^3','pi*G']
rel=mp.pslq(mp.matrix([I]+basis),tol=mp.mpf('1e-70'),maxcoeff=100000)
print('pslq:',rel)
if rel is not None and rel[0]!=0:
    print('  => I = ', -sum(mp.mpf(rel[k+1])*basis[k] for k in range(len(basis)))/rel[0])
    print('  z3 coeff =', mp.nstr(-mp.mpf(rel[1])/rel[0],10))
# 2) our c30*l0^3 for m=2
data=pickle.load(open('allcoeffs.pkl','rb'))
coeffs={k:sp.sympify(v) for k,v in data.items()}
_,_,syms=sym_G.build()
order=[syms[k] for k in ['s0','s1','c0','c1','tau1','tau2','tau3','gam1','gam2','gam3','S0','Sc','lS0']]
lam={k:sp.lambdify(order,coeffs[k],'mpmath') for k in coeffs}
def symvals(m,th):
    zeros=regular_zeros(m);sr=[mp.sin(th-b) for b in zeros]
    sv=[abs(x) for x in sr];cv=[mp.cos(th-b)*mp.sign(x) for b,x in zip(zeros,sr)]
    lv=[mp.log(s) for s in sv];N=range(1,m-1)
    tau=[mp.fsum([sv[i]*lv[i]**k for i in N]) for k in (1,2,3)]
    gam=[mp.fsum([cv[i]*(lv[i]**k+k*lv[i]**(k-1)) for i in N]) for k in (1,2,3)]
    S0=mp.fsum(sv);Sc=mp.fsum(cv)
    return [sv[0],sv[m-1],cv[0],cv[m-1],tau[0],tau[1],tau[2],gam[0],gam[1],gam[2],S0,Sc,mp.log(S0)],lv[0],lv[1]
for m in [2]:
    A=mp.pi/m
    J=mp.quad(lambda th: lam[(3,0)](*symvals(m,th)[0])*symvals(m,th)[1]**3,[0,A],method='tanh-sinh')
    print('m=2 int c30*l0^3 =',mp.nstr(J,30))
    rel=mp.pslq(mp.matrix([J]+basis),tol=mp.mpf('1e-68'),maxcoeff=100000)
    print('   pslq:',rel)
    if rel is not None and rel[0]!=0: print('   z3 coeff=',mp.nstr(-mp.mpf(rel[1])/rel[0],10))
