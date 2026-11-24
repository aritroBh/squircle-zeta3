import os, subprocess
d17 = "research/08_theory/D17_zeta5_polygons"
venv_py = "./.venv/bin/python"

script = r'''
import mpmath as mp, sys, time, json
mp.mp.dps = 60
print("dps:", mp.mp.dps)
from math import factorial
def regular_zeros(m):
    A=mp.pi/m; return [mp.mpf(0)]+[-(i*A) for i in range(1,m)]
def D5_pert(m, sx, sy, th, zeros):
    sraw=[mp.sin(th-b) for b in zeros]; sv=[abs(x) for x in sraw]
    cv=[mp.cos(th-b)*mp.sign(x) for b,x in zip(zeros,sraw)]
    ls=[mp.log(s) for s in sv]; ls[0]+=sx; ls[m-1]+=sy
    w=[1]*m
    T=[mp.fsum([w[i]*sv[i]*ls[i]**k for i in range(m)]) for k in range(6)]
    Tp=[mp.fsum([w[i]*cv[i] for i in range(m)])]+[mp.fsum([w[i]*cv[i]*(ls[i]**k+k*ls[i]**(k-1)) for i in range(m)]) for k in range(1,6)]
    S0=T[0]; a=[T[k]/(factorial(k)*S0) for k in range(6)]
    U=[mp.mpf(0)]*6
    for k in range(1,6): U[k]=a[k]-mp.fsum([mp.mpf(k-j)/k*a[j]*U[k-j] for j in range(1,k)])
    d=[T[k]/factorial(k) for k in range(6)]; n=[Tp[k]/factorial(k) for k in range(6)]
    Up=[mp.mpf(0)]*6; Up[0]=n[0]/d[0]
    for k in range(1,6): Up[k]=(n[k]-mp.fsum([d[j]*Up[k-j] for j in range(1,k+1)]))/d[0]
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
    for k in range(1,6): E[k]=(W[k]/B-mp.fsum([E[i]*E[k-i] for i in range(1,k)]))/2
    D0=mp.sqrt(B)/S0
    return D0*mp.fsum([rho[i]*E[5-i] for i in range(6)])
PTS=[(i,j) for i in range(6) for j in range(6) if i+j<=5]
TERMS=[(i,j) for i in range(6) for j in range(6) if i+j<=5]
MX=mp.matrix([[mp.mpf(sx)**i*mp.mpf(sy)**j for (i,j) in TERMS] for sx,sy in PTS])
LU = mp.LU_decomp(MX)
IDX = {t:k for k,t in enumerate(TERMS)}
i41,i32,i23,i14 = IDX[(4,1)],IDX[(3,2)],IDX[(2,3)],IDX[(1,4)]
def Gsum_pert(m, th, zeros):
    vals=mp.matrix([D5_pert(m, mp.mpf(sx), mp.mpf(sy), th, zeros) for sx,sy in PTS])
    cf=mp.lu_solve(LU, vals)
    return cf[i41]+cf[i32]+cf[i23]+cf[i14]
def C5_closed(m):
    h=mp.pi/(2*m); Sh=mp.sin(h); Ch=mp.cos(h)
    return 2880*m*Sh**4*Ch*(65*Sh**6-63*Sh**4+17*Sh**2-1)
out={}
t0=time.time()
for m in range(2,9):
    A=mp.pi/m; zeros=regular_zeros(m)
    J5 = mp.quad(lambda th: Gsum_pert(m,th,zeros), [0,A], method='tanh-sinh')
    C5_pert = 5760*m*J5
    C5_cl = C5_closed(m)
    res = abs(C5_pert-C5_cl)
    out[m]=dict(J5=mp.nstr(J5,50),C5_pert=mp.nstr(C5_pert,50),C5_closed=mp.nstr(C5_cl,50),resid=mp.nstr(res,3))
    print("m=%d (%.0fs): C5_pert=%s  C5_closed=%s  resid=%s"%(m,time.time()-t0,mp.nstr(C5_pert,20),mp.nstr(C5_cl,20),mp.nstr(res,2)))
json.dump(out, open("pert_extract_m2_8.json","w"))
print("saved pert_extract_m2_8.json")
'''.replace("D17DIR", d17)
with open(os.path.join(d17,"pert_extract.py"),"w") as f: f.write(script)
r = subprocess.run([venv_py, os.path.join(d17,"pert_extract.py")], capture_output=True, text=True, cwd=d17, timeout=1500)
print(r.stdout); print("ERR:", r.stderr[-2000:])
