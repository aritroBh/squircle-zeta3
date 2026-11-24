import json, time
from mpmath import mp, mpf, log, pi, sqrt, exp, factorial, sin, cos, matrix, lu_solve
t0=time.time()
mp.dps = 100
z0 = (-log(2)+pi*1j)/2
R = abs(z0); th = mp.arg(z0)
j120 = json.load(open("../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json"))
b_tau = [mpf(j120['eta'][k])/factorial(k) for k in range(121)]
b_d   = [mpf(j120['d'][k])/factorial(k) for k in range(121)]

def linfit(b, beta, k0, k1=120, ncorr=5, Rf=None, thf=None):
    Rf = R if Rf is None else Rf; thf = th if thf is None else thf
    rows=[]; y=[]
    for k in range(k0,k1+1):
        c = cos((k-mpf('0.5'))*thf); s = sin((k-mpf('0.5'))*thf)
        row=[c, s]
        for jj in range(1,ncorr+1): row += [c/k**jj, s/k**jj]
        rows.append(row); y.append(b[k]*Rf**k*k**(-beta))
    n = len(rows[0]); A = matrix(rows); Y = matrix(y)
    ATA = matrix(n,n); ATY = matrix(n,1)
    for i in range(n):
        for jj in range(n): ATA[i,jj] = sum(A[r,i]*A[r,jj] for r in range(len(rows)))
        ATY[i] = sum(A[r,i]*Y[r] for r in range(len(rows)))
    x = lu_solve(ATA, ATY)
    res = sqrt(sum((sum(A[r,i]*x[i] for i in range(n))-Y[r])**2 for r in range(len(rows)))/len(rows))
    return [x[i] for i in range(n)], res

print("== refined C fits, 5 corrections ==")
print("theory: tau C1 =", mp.nstr(-sqrt(3*pi*R/8),16), "  d C2 =", mp.nstr(sqrt(3*pi*R/4),16))
for k0 in [70,80,90,100,110]:
    xt, rt = linfit(b_tau, mpf('-1.5'), k0)
    xd, rd = linfit(b_d, mpf('-1.5'), k0)
    print(f"k0={k0}: tau C1={mp.nstr(xt[0],16)} (err {mp.nstr(abs(xt[0]+sqrt(3*pi*R/8)),3)}) | d C2={mp.nstr(xd[1],16)} (err {mp.nstr(abs(xd[1]-sqrt(3*pi*R/4)),3)})")

print("\n== 1D profiles for R and theta (5 corr, window 70-120) ==")
def profile(b, tag):
    print(f" {tag}:")
    for iR in [-30,-20,-10,-5,0,5,10,20,30]:
        Rf = R*(1+mpf(iR)/1000000)
        x, res = linfit(b, mpf('-1.5'), 70, 120, 5, Rf, None)
        print(f"   R*(1{iR:+d}e-6): resid={mp.nstr(res,4)}")
    for iT in [-30,-20,-10,-5,0,5,10,20,30]:
        thf = th*(1+mpf(iT)/1000000)
        x, res = linfit(b, mpf('-1.5'), 70, 120, 5, None, thf)
        print(f"   th*(1{iT:+d}e-6): resid={mp.nstr(res,4)}")
profile(b_tau, "tau")
