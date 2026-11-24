import json, time
from mpmath import mp, mpf, log, pi, sqrt, exp, factorial, sin, cos, matrix, lu_solve
t0=time.time()
mp.dps = 100
z0 = (-log(2)+pi*1j)/2
R = abs(z0); th = mp.arg(z0)
j120 = json.load(open("../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json"))
b_tau = [mpf(j120['eta'][k])/factorial(k) for k in range(121)]
b_d   = [mpf(j120['d'][k])/factorial(k) for k in range(121)]

def linfit(b, beta, k0, k1=120, ncorr=3, Rf=None, thf=None):
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

print("== d-series (theory: C1=0, C2=+sqrt(3 pi R/4) = %s)" % mp.nstr(sqrt(3*pi*R/4),14))
for k0 in [50,70,90,100]:
    x, res = linfit(b_d, mpf('-1.5'), k0)
    print(f"k0={k0}: C1={mp.nstr(x[0],6)} C2={mp.nstr(x[1],14)} resid={mp.nstr(res,3)}")

print("\n== (R,theta) scan, beta=-1.5 fixed, window 60-120 (log10 resid)")
def scan(b, tag):
    best=None
    for iR in range(-40,41,4):
        for iT in range(-40,41,4):
            Rf = R*(1+mpf(iR)/100000); thf = th*(1+mpf(iT)/100000)
            x, res = linfit(b, mpf('-1.5'), 60, 120, 3, Rf, thf)
            if best is None or res<best[0]: best=(res,Rf,thf)
    res,Rf,thf = best
    print(f"{tag}: best R={mp.nstr(Rf,12)} (rel err {mp.nstr(abs(Rf-R)/R,3)}), theta={mp.nstr(thf,12)} (rel err {mp.nstr(abs(thf-th)/th,3)}), resid={mp.nstr(res,3)}")
scan(b_tau, "tau-series")
scan(b_d, "d-series")
print(f"theory R={mp.nstr(R,12)}, theta={mp.nstr(th,12)}; A=1/R={mp.nstr(1/R,14)}")
