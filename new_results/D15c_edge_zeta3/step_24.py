
explore = r'''
import mpmath as mp
mp.mp.dps = 100
import sympy as sp
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
print("dps:",mp.mp.dps)
z3=mp.zeta(3)
for m in [2,3,4]:
    A=mp.pi/m
    coeffs=make_coeffs(regular_zeros(m),[1]*m,[0]*m)
    fullint=mp.quad(lambda th:coeffs(th)[3],[0,A],method='tanh-sinh')
    def cross(th):
        v=symvals(m,th); s0,s1=v[0],v[1]
        l0=mp.log(s0); l1=mp.log(s1)
        return f21(*v)*l0**2*l1 + f12(*v)*l0*l1**2
    crossint=mp.quad(cross,[0,A],method='tanh-sinh')
    crossG=mp.quad(lambda th: f21(*symvals(m,th))+f12(*symvals(m,th)),[0,A],method='tanh-sinh')
    Rint=fullint-crossint
    print(f"m={m}:")
    print(f"  fullint={mp.nstr(fullint,25)}")
    print(f"  [z3]fullint should be 2*z3*crossG; 2*crossG={mp.nstr(2*crossG,20)}")
    print(f"  crossint (numeric)={mp.nstr(crossint,25)}")
    print(f"  Rint=fullint-crossint={mp.nstr(Rint,25)}")
    print(f"  Rint/z3={mp.nstr(Rint/z3,20)}")
    # test: is fullint's z3-coeff = 2*crossG? i.e. does (fullint - 2*z3*crossG) have no z3?
    rem = fullint - 2*z3*crossG
    print(f"  fullint - 2*z3*crossG = {mp.nstr(rem,25)}")
'''
with open('step_12_explore.py','w') as f: f.write(explore)
import subprocess
r=subprocess.run(['./.venv/bin/python','step_12_explore.py'],capture_output=True,text=True,timeout=600)
print(r.stdout); print(r.stderr[-1500:] if r.stderr else "")
