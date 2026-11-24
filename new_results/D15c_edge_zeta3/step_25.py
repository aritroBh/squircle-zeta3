
decomp = r'''
import mpmath as mp
mp.mp.dps = 80
import sympy as sp, pickle
from poly_coeffs import regular_zeros
data = pickle.load(open('allcoeffs.pkl','rb'))
coeffs = {k: sp.sympify(v) for k,v in data.items()}
import sym_G
_,_,syms = sym_G.build()
order = [syms[k] for k in ['s0','s1','c0','c1','tau1','tau2','tau3','gam1','gam2','gam3','S0','Sc','lS0']]
lam = {k: sp.lambdify(order, coeffs[k], 'mpmath') for k in coeffs}
def symvals(m,th):
    zeros=regular_zeros(m); sr=[mp.sin(th-b) for b in zeros]
    sv=[abs(x) for x in sr]; cv=[mp.cos(th-b)*mp.sign(x) for b,x in zip(zeros,sr)]
    lv=[mp.log(s) for s in sv]; N=range(1,m-1)
    s0=sv[0];c0=cv[0];s1=sv[m-1];c1=cv[m-1]
    tau=[mp.fsum([sv[i]*lv[i]**k for i in N]) for k in (1,2,3)]
    gam=[mp.fsum([cv[i]*(lv[i]**k+k*lv[i]**(k-1)) for i in N]) for k in (1,2,3)]
    S0=mp.fsum(sv);Sc=mp.fsum(cv)
    return [s0,s1,c0,c1,tau[0],tau[1],tau[2],gam[0],gam[1],gam[2],S0,Sc,mp.log(S0)], lv[0], lv[1]
print("dps:",mp.mp.dps)
z3=mp.zeta(3); pi=mp.pi
# weight-3 basis + lower for m=2 (A=pi/2)
basis = [z3, pi**3, pi**2*mp.log(2), pi*mp.log(2)**2, mp.log(2)**3, pi*mp.catalan]
bn = ['z3','pi^3','pi^2log2','pi log2^2','log2^3','pi*G']
for m in [2]:
    A=mp.pi/m
    print(f"=== m={m}, component integrals I_ij = int c_ij l0^i l1^j, zeta(3) coeff via PSLQ ===")
    tot=mp.mpf(0)
    for (i,j) in [(3,0),(0,3),(2,1),(1,2),(2,0),(0,2),(1,1),(1,0),(0,1),(0,0)]:
        def f(th):
            v,l0,l1 = symvals(m,th)
            return lam[(i,j)](*v) * l0**i * l1**j
        I = mp.quad(f, [0,A], method='tanh-sinh')
        tot += I
        # PSLQ [I, basis]
        v = [I]+basis
        rel = mp.pslq(mp.matrix(v), tol=mp.mpf('1e-70'), maxcoeff=10000)
        zc = None
        if rel is not None:
            rel=list(rel)
            if rel[0]!=0:
                # normalize so coeff of I is -1: I = sum rel[1:]/(-rel[0]) * basis
                zc = -rel[1]/rel[0]
        print(f"  ({i},{j}): I={mp.nstr(I,18)}  z3coeff={mp.nstr(zc,14) if zc is not None else 'NOREL'}")
    print(f"  sum of components = {mp.nstr(tot,18)}")
'''
with open('step_13_decomp.py','w') as f: f.write(decomp)
import subprocess
r=subprocess.run(['./.venv/bin/python','step_13_decomp.py'],capture_output=True,text=True,timeout=900)
print(r.stdout); print(r.stderr[-1500:] if r.stderr else "")
