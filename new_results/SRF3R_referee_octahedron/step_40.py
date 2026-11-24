import json, time
from mpmath import mp, mpf, log, pi, sqrt, exp, factorial, sin, cos, matrix, lu_solve
t0=time.time()
mp.dps = 100
z0 = (-log(2)+pi*1j)/2
R = abs(z0); th = mp.arg(z0)
j120 = json.load(open("../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json"))
b_tau = [mpf(j120['eta'][k])/factorial(k) for k in range(121)]

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

print("R profile (tau, k0=70, 5 corr):")
for iR in [-30,-20,-10,-5,0,5,10,20,30]:
    Rf = R*(1+mpf(iR)/1000000)
    x, res = linfit(b_tau, mpf('-1.5'), 70, 120, 5, Rf, None)
    print(f"  R*(1{iR:+d}e-6): resid={mp.nstr(res,4)}")
print("theta profile:")
for iT in [-30,-20,-10,-5,0,5,10,20,30]:
    thf = th*(1+mpf(iT)/1000000)
    x, res = linfit(b_tau, mpf('-1.5'), 70, 120, 5, None, thf)
    print(f"  th*(1{iT:+d}e-6): resid={mp.nstr(res,4)}")
# beta profile with 5 corrections
print("beta profile (5 corr, k0=70):")
for ib in range(-1510,-1489,2):
    x, res = linfit(b_tau, mpf(ib)/1000, 70, 120, 5)
    print(f"  beta={ib/1000:.3f}: resid={mp.nstr(res,4)}")
