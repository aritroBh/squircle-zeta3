import mpmath as mp
mp.mp.dps = 100
import sympy as sp, json
from poly_coeffs import make_coeffs, regular_zeros
import sym_G
G21s,G12s,syms = sym_G.build()
order = [syms[k] for k in ['s0','s1','c0','c1','tau1','tau2','tau3','gam1','gam2','gam3','S0','Sc','lS0']]
f21 = sp.lambdify(order,G21s,'mpmath'); f12 = sp.lambdify(order,G12s,'mpmath')
def symvals(m,th):
    zeros=regular_zeros(m); sr=[mp.sin(th-b) for b in zeros]
    sv=[abs(x) for x in sr]; cv=[mp.cos(th-b)*mp.sign(x) for b,x in zip(zeros,sr)]
    lv=[mp.log(s) for s in sv]; N=range(1,m-1)
    s0=sv[0];c0=cv[0];s1=sv[m-1];c1=cv[m-1]
    tau=[mp.fsum([sv[i]*lv[i]**k for i in N]) for k in (1,2,3)]
    gam=[mp.fsum([cv[i]*(lv[i]**k+k*lv[i]**(k-1)) for i in N]) for k in (1,2,3)]
    S0=mp.fsum(sv);Sc=mp.fsum(cv)
    return [s0,s1,c0,c1,tau[0],tau[1],tau[2],gam[0],gam[1],gam[2],S0,Sc,mp.log(S0)]
print('dps:',mp.mp.dps)
z3=mp.zeta(3)
vals=json.load(open('vals_m2_8_200.json'))
def Cm(m):
    h=mp.pi/(2*m); s=mp.sin(h); return 48*m*s**4*mp.cos(h)*(3*s**2-1)
out={}
for m in range(2,9):
    A=mp.pi/m
    coeffs=make_coeffs(regular_zeros(m),[1]*m,[0]*m)
    fullint=mp.quad(lambda th:coeffs(th)[3],[0,A],method='tanh-sinh')
    crossG=mp.quad(lambda th: f21(*symvals(m,th))+f12(*symvals(m,th)),[0,A],method='tanh-sinh')
    P3_full = 12*m*fullint
    P3v = mp.mpf(vals[str(m)][3])
    rem_m = fullint - 2*z3*crossG
    out[m]=dict(fullint=mp.nstr(fullint,90), crossG=mp.nstr(crossG,90), rem=mp.nstr(rem_m,90))
    print('m=%d: crossG-C_m/(24m) resid=%s'%(m, mp.nstr(abs(crossG-Cm(m)/(24*m)),3)))
    print('    12m*fullint - vals P3 resid=%s'%mp.nstr(abs(P3_full-P3v),3))
    print('    rem_m=%s'%mp.nstr(rem_m,28))
json.dump(out,open('rem_m.json','w'))
