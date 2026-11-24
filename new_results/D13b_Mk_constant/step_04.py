
import mpmath as mp
mp.mp.dps = 100

def load(fn):
    out = {}
    for line in open(fn):
        parts = line.split("\t")
        k = int(parts[0]); mid = parts[1].strip()
        out[k] = mp.mpf(mid.split(" ")[0].strip("[]"))
    return out

d = load("mk_deg9.tsv")
spi = mp.sqrt(mp.pi)
R = {k: d[k]*k*mp.exp(mp.sqrt(mp.pi*k))/spi for k in range(2,401)}
psi = {k: mp.sqrt(mp.pi*k) for k in range(2,401)}

def fit(kmin, J, kmax=400):
    # model: R_k = Re[ e^{i psi_k} * sum_{j=0}^{J} c_j k^{-j/2} ], c_j = a_j + i b_j
    ks = list(range(kmin, kmax+1))
    n = len(ks); m = 2*(J+1)
    A = mp.matrix(n, m); y = mp.matrix(n,1)
    for i,k in enumerate(ks):
        for j in range(J+1):
            x = k**(-mp.mpf(j)/2)
            A[i, 2*j]   = mp.cos(psi[k])*x
            A[i, 2*j+1] = -mp.sin(psi[k])*x
        y[i] = R[k]
    # normal equations
    At = A.T
    G = At*A; rhs = At*y
    sol = mp.lu_solve(G, rhs)
    cs = [(sol[2*j], sol[2*j+1]) for j in range(J+1)]
    # residual rms
    rms = mp.mpf(0)
    for i,k in enumerate(ks):
        pred = sum((A[i,2*j]*sol[2*j] + A[i,2*j+1]*sol[2*j+1]) for j in range(J+1))
        rms += (pred-y[i])**2
    rms = mp.sqrt(rms/n)
    return cs, rms

def report(kmin, J):
    cs, rms = fit(kmin, J)
    a0,b0 = cs[0]
    A0 = mp.hypot(a0,b0); phi0 = mp.atan2(b0,a0)
    if phi0 < 0: phi0 += 2*mp.pi
    # delta method for A1, phi1
    a1,b1 = cs[1]
    A1 = (a0*a1+b0*b1)/A0
    phi1 = (a0*b1-b0*a1)/A0**2
    a2,b2 = cs[2]
    A2 = (a0*a2+b0*b2)/A0
    phi2 = (a0*b2-b0*a2)/A0**2
    print(f"kmin={kmin:4d} J={J}  A0={mp.nstr(A0,18)}  phi0={mp.nstr(phi0,18)} ({mp.nstr(phi0/mp.pi,10)} pi)")
    print(f"             A1={mp.nstr(A1,12)} phi1={mp.nstr(phi1,12)}   A2={mp.nstr(A2,10)} phi2={mp.nstr(phi2,10)}  rms={mp.nstr(rms,3)}")

for J in [2,3,4]:
    for kmin in [60, 100, 150, 200]:
        report(kmin, J)
    print()
