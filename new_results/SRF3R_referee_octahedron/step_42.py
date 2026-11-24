import json
from mpmath import mp, mpf, log, pi, sqrt, zeta, factorial, cos, sin, matrix, lu_solve
mp.dps = 80
z0 = (-log(2)+pi*1j)/2; R = abs(z0); th = mp.arg(z0)

# exact low-order d-coefficients from paper thm:oct
j120 = json.load(open("../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json"))
d = [mpf(x) for x in j120['d']]
S1 = 4*sqrt(3)
ex2 = 1-pi**2/9
ex3 = 1+pi**2/9-16*zeta(3)/9
ex4 = 1-2*pi**2/9+32*zeta(3)/9-16*pi**4/405
print("d^2:", mp.nstr(d[2]-ex2,5), "| d^3:", mp.nstr(d[3]-ex3,5), "| d^4:", mp.nstr(d[4]-ex4,5), " (should be ~0)")

# w-series: log[(2-p)S/(4 sqrt3)] = sum w_k tau^k ; paper values
wjson = json.load(open("../SRF3_octahedron_radius_proof/w_series_K120.json"))
print("w_series keys:", wjson.keys() if isinstance(wjson, dict) else "list len "+str(len(wjson)))
wlist = wjson['w'] if isinstance(wjson, dict) and 'w' in wjson else wjson
w = [mpf(x) for x in wlist]
wex = [None, None, -pi**2/9, -16*zeta(3)/9, -37*pi**4/810, -160*zeta(5)/27, 16*zeta(3)**2/81 - 572*pi**6/76545]
for k in range(2,7):
    print(f"w_{k}: diff vs paper exact = {mp.nstr(w[k]-wex[k],5)}")

# w-series leading law: fit C_w, compare with C_tau (should be identical, claim C_w/C_tau=1)
bw = [w[k]/factorial(k) for k in range(len(w))]
def linfit(b, beta, k0, k1=120, ncorr=5):
    rows=[]; y=[]
    for k in range(k0,k1+1):
        c = cos((k-mpf('0.5'))*th); s = sin((k-mpf('0.5'))*th)
        row=[c, s]
        for jj in range(1,ncorr+1): row += [c/k**jj, s/k**jj]
        rows.append(row); y.append(b[k]*R**k*k**(-beta))
    n = len(rows[0]); A = matrix(rows); Y = matrix(y)
    ATA = matrix(n,n); ATY = matrix(n,1)
    for i in range(n):
        for jj in range(n): ATA[i,jj] = sum(A[r,i]*A[r,jj] for r in range(len(rows)))
        ATY[i] = sum(A[r,i]*Y[r] for r in range(len(rows)))
    x = lu_solve(ATA, ATY)
    return [x[i] for i in range(n)]
xw = linfit(bw, mpf('-1.5'), 80)
print("w-series: C1 =", mp.nstr(xw[0],16), " C2 =", mp.nstr(xw[1],6))
print("tau C1 theory:", mp.nstr(-sqrt(3*pi*R/8),16), " ratio:", mp.nstr(xw[0]/(-sqrt(3*pi*R/8)),12))
