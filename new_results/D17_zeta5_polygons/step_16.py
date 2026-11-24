import os
d17 = "research/08_theory/D17_zeta5_polygons"

builder = r'''
import sympy as sp

def build_D(kmax=5, reduced=False):
    """Build D_0..D_kmax symbolically via explicit recursions (fast poly arithmetic).
    reduced=True sets tau_k=gam_k=lS0=0 (only edge symbols s0,s1,c0,c1,S0,Sc,L0,L1)."""
    s0,s1,c0,c1 = sp.symbols('s0 s1 c0 c1')
    L0,L1 = sp.symbols('L0 L1')
    S0 = sp.Symbol('S0'); Sc = sp.Symbol('Sc'); lS0 = sp.Symbol('lS0')
    if reduced:
        tau = [sp.Integer(0)]*kmax
        gam = [sp.Integer(0)]*kmax
        lS0v = sp.Integer(0)
    else:
        tau = list(sp.symbols('tau1:'+str(kmax+1)))
        gam = list(sp.symbols('gam1:'+str(kmax+1)))
        lS0v = lS0
    # Moments
    T = [S0]; Tp = [Sc]
    for k in range(1,kmax+1):
        T.append(tau[k-1] + s0*L0**k + s1*L1**k)
        Tp.append(gam[k-1] + c0*(L0**k + k*L0**(k-1)) + c1*(L1**k + k*L1**(k-1)))
    # a_k = T_k/(k! S0)
    import math
    a = [sp.expand(T[k]/(math.factorial(k)*S0)) for k in range(kmax+1)]  # a_0 = 1
    # U_k: log(1+ (S/S0-1)); U_k = a_k - (1/k) sum_{j=1}^{k-1} (k-j) a_j U_{k-j}
    U = [None]*(kmax+1)
    for k in range(1,kmax+1):
        s = a[k]
        for j in range(1,k):
            s -= sp.Rational(k-j,k)*a[j]*U[k-j]
        U[k] = sp.expand(s)
    # Up_k: d/dth log S coefficients; d_k = T_k/k!, n_k = Tp_k/k!
    d = [T[k]/math.factorial(k) for k in range(kmax+1)]
    n = [Tp[k]/math.factorial(k) for k in range(kmax+1)]
    Up = [None]*(kmax+1)
    Up[0] = sp.expand(n[0]/d[0])  # = Sc/S0 = L
    for k in range(1,kmax+1):
        s = n[k]
        for j in range(1,k+1):
            s -= d[j]*Up[k-j]
        Up[k] = sp.expand(s/d[0])
    # R_k = -sum_{j=0}^k (-1)^j u_{k-j}; u_0=lS0
    u = [lS0v]+U[1:]
    R = [None]*(kmax+1)
    for k in range(1,kmax+1):
        R[k] = sp.expand(-sum((-1)**j*u[k-j] for j in range(k+1)))
    # rho_k = [e^k] exp(sum R_j e^j); rho_0=1; rho_k=(1/k) sum_{i=1}^k i R_i rho_{k-i}
    rho = [sp.Integer(1)]*(kmax+1)
    for k in range(1,kmax+1):
        s = sum(i*R[i]*rho[k-i] for i in range(1,k+1))
        rho[k] = sp.expand(sp.Rational(1,k)*s)
    # Q_k = -sum_{j=0}^k (-1)^j Up_{k-j}; Up_0 = L
    Q = [None]*(kmax+1)
    for k in range(kmax+1):
        Q[k] = sp.expand(-sum((-1)**j*Up[k-j] for j in range(k+1)))
    # B = 1+Q0^2 ; W_k = sum_{i=0}^k Q_i Q_{k-i}
    B = sp.expand(1+Q[0]**2)
    W = [None]*(kmax+1)
    for k in range(1,kmax+1):
        W[k] = sp.expand(sum(Q[i]*Q[k-i] for i in range(k+1)))
    # E_k: F=sqrt(1+(W/B)); 2 E_k = W_k/B - sum_{i=1}^{k-1} E_i E_{k-i}
    E = [sp.Integer(1)]*(kmax+1)
    for k in range(1,kmax+1):
        s = W[k]/B
        for i in range(1,k):
            s -= E[i]*E[k-i]
        E[k] = sp.expand(s/2)
    D0 = sp.sqrt(B)/S0
    D = [None]*(kmax+1)
    for k in range(kmax+1):
        D[k] = sp.expand(D0*sum(rho[i]*E[k-i] for i in range(k+1)))
    syms = dict(s0=s0,s1=s1,c0=c0,c1=c1,L0=L0,L1=L1,S0=S0,Sc=Sc,lS0=lS0,B=B,D0=D0)
    return D, syms

if __name__=='__main__':
        import time, sys
        km = int(sys.argv[1]) if len(sys.argv)>1 else 3
        red = (len(sys.argv)>2 and sys.argv[2]=='red')
        t0=time.time()
        D,syms = build_D(km, reduced=red)
        print(f"built D0..D{km} reduced={red} in {time.time()-t0:.1f}s")
        for k in range(km+1):
            print("D%d terms:"%k, len(sp.expand(D[k]).as_ordered_terms()))
'''
with open(os.path.join(d17, "sym_D5.py"), "w") as f:
    f.write(builder)
print("saved sym_D5.py (recursion-based)")

import subprocess
venv_py = "./.venv/bin/python"
# First sanity: order 3 reduced, and compare structure to sym_G
r = subprocess.run([venv_py, os.path.join(d17,"sym_D5.py"), "3", "red"], capture_output=True, text=True, cwd=d17, timeout=300)
print(r.stdout); print("ERR:", r.stderr[-1500:])
