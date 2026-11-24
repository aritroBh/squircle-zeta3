import json
from mpmath import mp, mpf, log, pi, sqrt, exp, factorial, sin, cos, matrix, lu_solve
mp.dps = 100
z0 = (-log(2)+pi*1j)/2
R = abs(z0); th = mp.arg(z0)
j120 = json.load(open("../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json"))
b_tau = [mpf(j120['eta'][k])/factorial(k) for k in range(121)]
b_d   = [mpf(j120['d'][k])/factorial(k) for k in range(121)]

def linfit(b, beta, k0, k1=120, ncorr=5, extra_angle=None):
    rows=[]; y=[]
    for k in range(k0,k1+1):
        c = cos((k-mpf('0.5'))*th); s = sin((k-mpf('0.5'))*th)
        row=[c, s]
        for jj in range(1,ncorr+1): row += [c/k**jj, s/k**jj]
        if extra_angle is not None:
            cv = cos((k-mpf('0.5'))*extra_angle); sv = sin((k-mpf('0.5'))*extra_angle)
            row += [cv, sv]
        rows.append(row); y.append(b[k]*R**k*k**(-beta))
    n = len(rows[0]); A = matrix(rows); Y = matrix(y)
    ATA = matrix(n,n); ATY = matrix(n,1)
    for i in range(n):
        for jj in range(n): ATA[i,jj] = sum(A[r,i]*A[r,jj] for r in range(len(rows)))
        ATY[i] = sum(A[r,i]*Y[r] for r in range(len(rows)))
    x = lu_solve(ATA, ATY)
    return [x[i] for i in range(n)]

# residual table: model vs coefficients, k=60..120 (relative), using theory constants only
C1t = -sqrt(3*pi*R/8); C2d = sqrt(3*pi*R/4)
worst_t = 0; worst_d = 0
print(" k | rel resid tau | rel resid d (model = pure leading law, no corrections)")
for k in range(60,121,10):
    mt = C1t*R**(-k)*k**mpf('-1.5')*cos((k-mpf('0.5'))*th)
    md = C2d*R**(-k)*k**mpf('-1.5')*sin((k-mpf('0.5'))*th)
    rt_ = abs(mt-b_tau[k])/abs(b_tau[k]); rd_ = abs(md-b_d[k])/abs(b_d[k])
    worst_t = max(worst_t, rt_); worst_d = max(worst_d, rd_)
    print(f"{k:4d} | {mp.nstr(rt_,3)} | {mp.nstr(rd_,3)}")
print("worst in table:", mp.nstr(worst_t,3), mp.nstr(worst_d,3))

# vertex angle pi - theta = 1.3536: add as free component, window 50-120
vth = pi - th
x = linfit(b_tau, mpf('-1.5'), 50, 120, 5, extra_angle=vth)
amp = sqrt(x[-2]**2 + x[-1]**2)
print(f"\nvertex angle {mp.nstr(vth,8)}: fitted amplitude = {mp.nstr(amp,4)} vs edge |C1| = {mp.nstr(abs(x[0]),8)}")
print(f"amplitude ratio = {mp.nstr(amp/abs(x[0]),4)} (SRF3 claim: < 2e-3)")
x2 = linfit(b_d, mpf('-1.5'), 50, 120, 5, extra_angle=vth)
amp2 = sqrt(x2[-2]**2 + x2[-1]**2)
print(f"d-series vertex amplitude ratio = {mp.nstr(amp2/abs(x2[1]),4)}")
