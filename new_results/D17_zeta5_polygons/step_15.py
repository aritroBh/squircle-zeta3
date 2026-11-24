import os, json
d17 = "research/08_theory/D17_zeta5_polygons"
os.makedirs(d17, exist_ok=True)

# Save universal weights derived + verified
weights = {
  "weight3": {"(2,1)": {"z3": 2}, "(1,2)": {"z3": 2}},
  "weight5": {
    "(4,1)": {"z5": 24, "z2z3": 0},
    "(1,4)": {"z5": 24, "z2z3": 0},
    "(3,2)": {"z5": 24, "z2z3": -12},
    "(2,3)": {"z5": 24, "z2z3": -12},
  },
  "rule": "[zeta(n)] int_0^1 log^a t log^b(1-t) dt = single-zeta weight; z5 weight = a!b!*24/(a!b!)=24 for all cross",
  "proof_sketch": "log B(1+u,1+v)=logG(1+u)+logG(1+v)-logG(2+u+v); cross u^a v^b coeff (a+b=5) = -psi4(2)/(a!b!)=24(z5-1)/(a!b!); moment=a!b!*[coeff in exp]; single z5 only from direct L term"
}
with open(os.path.join(d17, "universal_weights.json"), "w") as f:
    json.dump(weights, f, indent=2)
print("saved universal_weights.json")

# Now build the order-5 symbolic machinery, extending sym_G.py
builder = r'''
import sympy as sp

def build_D(kmax=5):
    # symbols
    s0,s1,c0,c1 = sp.symbols('s0 s1 c0 c1')
    L0,L1 = sp.symbols('L0 L1')
    tau = sp.symbols('tau1:'+str(kmax+1))  # tau1..tau5
    gam = sp.symbols('gam1:'+str(kmax+1))  # gam1..gam5
    S0 = sp.Symbol('S0'); Sc = sp.Symbol('Sc'); lS0 = sp.Symbol('lS0')
    e = sp.Symbol('e')
    # Moments T_k and theta-derivatives Tp_k
    T = [S0]
    Tp = [Sc]
    for k in range(1,kmax+1):
        T.append(tau[k-1] + s0*L0**k + s1*L1**k)
        Tp.append(gam[k-1] + c0*(L0**k + k*L0**(k-1)) + c1*(L1**k + k*L1**(k-1)))
    # S(e) = sum_k T_k e^k/k!  ; log S = lS0 + sum U_k e^k
    Sser = sum(T[k]/sp.factorial(k)*e**k for k in range(kmax+1))
    logS = sp.log(Sser)  # will expand via series
    # Use series expansion of log
    logS_exp = sp.series(sp.log(Sser), e, 0, kmax+1).removeO().expand()
    U = [sp.expand(logS_exp).coeff(e,k) for k in range(1,kmax+1)]
    # S'/S : numerator Nser = sum Tp_k e^k/k!, denom Sser
    Nser = sum(Tp[k]/sp.factorial(k)*e**k for k in range(kmax+1))
    dlogS = sp.series(Nser/Sser, e, 0, kmax+1).removeO().expand()
    Up = [sp.expand(dlogS).coeff(e,k) for k in range(1,kmax+1)]
    # R_k = -sum_{j=0}^k (-1)^j u_{k-j}, u_0=lS0
    u = [lS0]+U
    R = [0]*(kmax+1)
    for k in range(1,kmax+1):
        R[k] = sp.expand(-sum((-1)**j*u[k-j] for j in range(k+1)))
    # rho_k = [e^k] exp(sum R_j e^j)
    Rser = sum(R[k]*e**k for k in range(1,kmax+1))
    expR = sp.series(sp.exp(Rser), e, 0, kmax+1).removeO().expand()
    rho = [1]+[sp.expand(expR).coeff(e,k) for k in range(1,kmax+1)]
    # Q_k = d/dth R_k = -sum_{j=0}^k (-1)^j Up_{k-j}, Up_0 = L = Sc/S0
    L = Sc/S0
    up0 = [L]+Up
    Q = [0]*(kmax+1)
    Q[0] = -L
    for k in range(1,kmax+1):
        Q[k] = sp.expand(-sum((-1)**j*up0[k-j] for j in range(k+1)))
    # W_k = sum_{i=0}^k Q_i Q_{k-i}
    B = 1+Q[0]**2
    W = [0]*(kmax+1)
    for k in range(1,kmax+1):
        W[k] = sp.expand(sum(Q[i]*Q[k-i] for i in range(k+1)))
    # E_k = [e^k] sqrt(1 + (sum W_j e^j)/B)
    y = sum(W[k]/B*e**k for k in range(1,kmax+1))
    sq = sp.series(sp.sqrt(1+y), e, 0, kmax+1).removeO().expand()
    E = [1]+[sp.expand(sq).coeff(e,k) for k in range(1,kmax+1)]
    D0 = sp.sqrt(B)/S0
    # D_k = D0 * sum_{i=0}^k rho_i E_{k-i}
    D = [0]*(kmax+1)
    for k in range(kmax+1):
        D[k] = sp.expand(D0*sum(rho[i]*E[k-i] for i in range(k+1)))
    syms = dict(s0=s0,s1=s1,c0=c1,c1=c1,L0=L0,L1=L1,S0=S0,Sc=Sc,lS0=lS0,e=e,
                tau=tau,gam=gam)
    return D, syms

if __name__=='__main__':
    import time
    t0=time.time()
    D,syms = build_D(5)
    print("built D0..D5 in", time.time()-t0,"s")
    for k in range(6):
        print("D%d terms:"%k, len(sp.expand(D[k]).as_ordered_terms()))
'''
with open(os.path.join(d17, "sym_D5.py"), "w") as f:
    f.write(builder)
print("saved sym_D5.py")

# Test the build
import subprocess, sys
venv_py = "./.venv/bin/python"
r = subprocess.run([venv_py, os.path.join(d17,"sym_D5.py")], capture_output=True, text=True, cwd=d17)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
