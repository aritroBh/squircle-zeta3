
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

# pointwise phasor estimates
def est(k):
    dk = psi[k+1]-psi[k]
    re = R[k]
    im = (R[k+1]-R[k]*mp.cos(dk))/mp.sin(dk)
    A = mp.hypot(re, im)
    th = -mp.atan2(im, re)
    ph = th - psi[k]
    ph = mp.fmod(ph, 2*mp.pi)
    if ph < 0: ph += 2*mp.pi
    return A, ph

def polyfit(vals, kmin, kmax, J):
    ks = [k for k in range(kmin, kmax+1)]
    n = len(ks); m = J+1
    A = mp.matrix(n,m); y = mp.matrix(n,1)
    for i,k in enumerate(ks):
        for j in range(m):
            A[i,j] = k**(-mp.mpf(j)/2)
        y[i] = vals[k]
    G = A.T*A; rhs = A.T*y
    return mp.lu_solve(G, rhs)

Avals = {}; pvals = {}
for k in range(50, 400):
    A, ph = est(k)
    Avals[k] = A; pvals[k] = ph

print("=== fits of pointwise estimates:  v_k = c0 + c1 k^{-1/2} + c2 k^{-1} + c3 k^{-3/2} (+c4 k^{-2}) ===")
for J in [2,3,4]:
    for kmin in [60, 100, 150, 200, 250]:
        ca = polyfit(Avals, kmin, 399, J)
        cp = polyfit(pvals, kmin, 399, J)
        print(f"J={J} kmin={kmin:3d}:  A0={mp.nstr(ca[0],16)}  phi0={mp.nstr(cp[0],16)} ({mp.nstr(cp[0]/mp.pi,12)} pi)   A1={mp.nstr(ca[1],10)} phi1={mp.nstr(cp[1],10)}")
    print()
