import os, subprocess
d17 = "research/08_theory/D17_zeta5_polygons"
venv_py = "./.venv/bin/python"

poly5 = r'''
import mpmath as mp
from math import factorial

def make_coeffs5(zeros, weights, lams):
    """Numerical D_0..D_5 via recursions (extends D15b poly_coeffs to order 5)."""
    def coeffs(th):
        sraw = [mp.sin(th - b) for b in zeros]
        sv = [abs(x) for x in sraw]
        cv = [mp.cos(th - b)*mp.sign(x) for b, x in zip(zeros, sraw)]
        ls = [mp.log(s) + lam for s, lam in zip(sv, lams)]
        w = weights
        T = [mp.fsum([w[i]*sv[i]*ls[i]**k for i in range(len(w))]) for k in range(6)]
        Tp = [mp.fsum([w[i]*cv[i] for i in range(len(w))])] + \
             [mp.fsum([w[i]*cv[i]*(ls[i]**k + k*ls[i]**(k-1)) for i in range(len(w))]) for k in range(1,6)]
        S0 = T[0]
        # a_k = T_k/(k! S0)
        a = [T[k]/(factorial(k)*S0) for k in range(6)]
        # U_k: log(S/S0); U_k = a_k - (1/k) sum_{j=1}^{k-1}(k-j) a_j U_{k-j}
        U = [mp.mpf(0)]*6
        for k in range(1,6):
            s = a[k]
            for j in range(1,k):
                s -= mp.mpf(k-j)/k*a[j]*U[k-j]
            U[k] = s
        # Up_k: d/dth log S; d_k=T_k/k!, n_k=Tp_k/k!
        d = [T[k]/factorial(k) for k in range(6)]
        n = [Tp[k]/factorial(k) for k in range(6)]
        Up = [mp.mpf(0)]*6
        Up[0] = n[0]/d[0]
        for k in range(1,6):
            s = n[k]
            for j in range(1,k+1):
                s -= d[j]*Up[k-j]
            Up[k] = s/d[0]
        # R_k = -sum_{j=0}^k (-1)^j u_{k-j}; u_0 = log S0
        u = [mp.log(S0)] + [U[k] for k in range(1,6)]
        R = [mp.mpf(0)]*6
        for k in range(1,6):
            R[k] = -mp.fsum([(-1)**j*u[k-j] for j in range(k+1)])
        # rho_k = (1/k) sum_{i=1}^k i R_i rho_{k-i}
        rho = [mp.mpf(1)]*6
        for k in range(1,6):
            rho[k] = mp.fsum([i*R[i]*rho[k-i] for i in range(1,k+1)])/k
        # Q_k = -sum_{j=0}^k (-1)^j Up_{k-j}
        Q = [mp.mpf(0)]*6
        for k in range(6):
            Q[k] = -mp.fsum([(-1)**j*Up[k-j] for j in range(k+1)])
        B = 1+Q[0]**2
        W = [mp.mpf(0)]*6
        for k in range(1,6):
            W[k] = mp.fsum([Q[i]*Q[k-i] for i in range(k+1)])
        E = [mp.mpf(1)]*6
        for k in range(1,6):
            s = W[k]/B
            for i in range(1,k):
                s -= E[i]*E[k-i]
            E[k] = s/2
        D0 = mp.sqrt(B)/S0
        return [D0*mp.fsum([rho[i]*E[k-i] for i in range(k+1)]) for k in range(6)]
    return coeffs

def regular_zeros(m):
    A = mp.pi/m
    return [mp.mpf(0)] + [-(i*A) for i in range(1,m)]

def regular_P_coeffs5(m, kmax=5):
    A = mp.pi/m
    coeffs = make_coeffs5(regular_zeros(m), [1]*m, [0]*m)
    x = A/2
    return [4*m*mp.quad(lambda th: coeffs(th)[k], [0, x], method='tanh-sinh') for k in range(kmax+1)]

def edge_arc_coeffs5(m, kmax=5):
    A = mp.pi/m
    coeffs = make_coeffs5(regular_zeros(m), [1]*m, [0]*m)
    return [2*m*mp.quad(lambda th: coeffs(th)[k], [0, A], method='tanh-sinh') for k in range(kmax+1)]
'''
with open(os.path.join(d17,"poly_coeffs5.py"),"w") as f: f.write(poly5)
print("saved poly_coeffs5.py")

# Validate D3 against existing poly_coeffs (order 3), m=2..5
val = r'''
import mpmath as mp, sys
mp.mp.dps = 80
sys.path.insert(0, "D17DIR")
sys.path.insert(0, "D15BDIR")
import poly_coeffs5 as p5
import poly_coeffs as p3
print("dps:", mp.mp.dps)
for m in [2,3,4,5]:
    A = mp.pi/m
    c5 = p5.make_coeffs5(p5.regular_zeros(m), [1]*m, [0]*m)
    c3 = p3.make_coeffs(p3.regular_zeros(m), [1]*m, [0]*m)
    th = mp.mpf('0.37')*A
    v5 = c5(th); v3 = c3(th)
    res = max(abs(v5[k]-v3[k]) for k in range(4))
    print("m=%d th=%s : max |D5rec-D3old| for k=0..3 = %s"%(m, mp.nstr(th,6), mp.nstr(res,3)))
'''.replace("D17DIR", d17).replace("D15BDIR", "research/08_theory/D15b_polygon_law")
with open(os.path.join(d17,"val_poly5.py"),"w") as f: f.write(val)
r = subprocess.run([venv_py, os.path.join(d17,"val_poly5.py")], capture_output=True, text=True, cwd=d17, timeout=300)
print(r.stdout); print("ERR:", r.stderr[-2000:])
