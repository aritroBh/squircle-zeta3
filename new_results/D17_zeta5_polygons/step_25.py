import os, subprocess
d17 = "research/08_theory/D17_zeta5_polygons"
venv_py = "./.venv/bin/python"

script = r'''
import sympy as sp, sys, time
import mpmath as mp
mp.mp.dps = 70
sys.path.insert(0, "D17DIR")
import sym_D5
from math import factorial

# ---- reduced symbolic cross G's ----
D, syms = sym_D5.build_D(5, reduced=True)
L0,L1 = syms['L0'], syms['L1']
s0,s1,c0,c1,S0,Sc = syms['s0'],syms['s1'],syms['c0'],syms['c1'],syms['S0'],syms['Sc']
D5 = sp.expand(D[5])
symG = {}
for (a,b) in [(4,1),(3,2),(2,3),(1,4)]:
    symG[(a,b)] = sp.lambdify((s0,s1,c0,c1,S0,Sc), sp.expand(D5.coeff(L0,a).coeff(L1,b)), 'mpmath')

# ---- perturbative extractor (full, all edges) ----
def regular_zeros(m):
    A=mp.pi/m; return [mp.mpf(0)]+[-(i*A) for i in range(1,m)]
def D5_pert(m, sx, sy, th):
    zeros = regular_zeros(m)
    sraw=[mp.sin(th-b) for b in zeros]; sv=[abs(x) for x in sraw]
    cv=[mp.cos(th-b)*mp.sign(x) for b,x in zip(zeros,sraw)]
    ls=[mp.log(s) for s in sv]; ls[0]+=sx; ls[m-1]+=sy
    w=[1]*m
    T=[mp.fsum([w[i]*sv[i]*ls[i]**k for i in range(m)]) for k in range(6)]
    Tp=[mp.fsum([w[i]*cv[i] for i in range(m)])]+[mp.fsum([w[i]*cv[i]*(ls[i]**k+k*ls[i]**(k-1)) for i in range(m)]) for k in range(1,6)]
    S0=T[0]; a=[T[k]/(factorial(k)*S0) for k in range(6)]
    U=[mp.mpf(0)]*6
    for k in range(1,6):
        U[k]=a[k]-mp.fsum([mp.mpf(k-j)/k*a[j]*U[k-j] for j in range(1,k)])
    d=[T[k]/factorial(k) for k in range(6)]; n=[Tp[k]/factorial(k) for k in range(6)]
    Up=[mp.mpf(0)]*6; Up[0]=n[0]/d[0]
    for k in range(1,6):
        Up[k]=(n[k]-mp.fsum([d[j]*Up[k-j] for j in range(1,k+1)]))/d[0]
    u=[mp.log(S0)]+[U[k] for k in range(1,6)]
    R=[mp.mpf(0)]*6
    for k in range(1,6): R[k]=-mp.fsum([(-1)**j*u[k-j] for j in range(k+1)])
    rho=[mp.mpf(1)]*6
    for k in range(1,6): rho[k]=mp.fsum([i*R[i]*rho[k-i] for i in range(1,k+1)])/k
    Q=[mp.mpf(0)]*6
    for k in range(6): Q[k]=-mp.fsum([(-1)**j*Up[k-j] for j in range(k+1)])
    B=1+Q[0]**2
    W=[mp.mpf(0)]*6
    for k in range(1,6): W[k]=mp.fsum([Q[i]*Q[k-i] for i in range(k+1)])
    E=[mp.mpf(1)]*6
    for k in range(1,6):
        E[k]=(W[k]/B-mp.fsum([E[i]*E[k-i] for i in range(1,k)]))/2
    D0=mp.sqrt(B)/S0
    return D0*mp.fsum([rho[i]*E[5-i] for i in range(6)])
PTS=[(i,j) for i in range(6) for j in range(6) if i+j<=5]
TERMS=[(i,j) for i in range(6) for j in range(6) if i+j<=5]
MX=mp.matrix([[mp.mpf(sx)**i*mp.mpf(sy)**j for (i,j) in TERMS] for sx,sy in PTS])
MXinv = MX**-1
def pertG(m, th):
    vals=mp.matrix([D5_pert(m, mp.mpf(sx), mp.mpf(sy), th) for sx,sy in PTS])
    cf=MXinv*vals
    dct=dict(zip(TERMS, list(cf)))
    return dct

def symvals(m, th):
    # s0,s1,c0,c1,S0,Sc for the edge (full sums over all m edges)
    zeros=regular_zeros(m); sr=[mp.sin(th-b) for b in zeros]
    sv=[abs(x) for x in sr]; cv=[mp.cos(th-b)*mp.sign(x) for b,x in zip(zeros,sr)]
    return (sv[0], sv[m-1], cv[0], cv[m-1], mp.fsum(sv), mp.fsum(cv))

print("dps:", mp.mp.dps)
for m in [3,4]:
    A=mp.pi/m
    for frac in ['0.31','0.57']:
        th=mp.mpf(frac)*A
        pg=pertG(m,th)
        svv=symvals(m,th)
        for (a,b) in [(4,1),(3,2),(2,3),(1,4)]:
            sg=symG[(a,b)](*svv)
            res=abs(pg[(a,b)]-sg)
            rel=res/max(mp.mpf(1),abs(sg))
            print("m=%d th=%sA G%d%d: pert=%s sym=%s relres=%s"%(m,frac,a,b,mp.nstr(pg[(a,b)],10),mp.nstr(sg,10),mp.nstr(rel,2)))
'''.replace("D17DIR", d17)
with open(os.path.join(d17,"check_indep.py"),"w") as f: f.write(script)
r = subprocess.run([venv_py, os.path.join(d17,"check_indep.py")], capture_output=True, text=True, cwd=d17, timeout=1200)
print(r.stdout); print("ERR:", r.stderr[-2000:])
