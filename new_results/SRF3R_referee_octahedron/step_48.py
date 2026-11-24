import json
from mpmath import mp, mpf, log, pi, sqrt, factorial, cos, sin, matrix, lu_solve
mp.dps = 80
z0 = (-log(2)+pi*1j)/2; R = abs(z0); th = mp.arg(z0)
j120 = json.load(open("../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json"))
A_tau = [mpf(x) for x in j120['eta']]
A_d   = [mpf(x) for x in j120['d']]
# mixing identity: b_k(d) = sum_{j=1}^k a_j(tau) (-1)^{k-j} C(k-1, j-1)
from math import comb
worst = mpf(0)
for k in range(1, 121):
    s = mpf(0)
    for jj in range(1, k+1):
        s += A_tau[jj]*((-1)**(k-jj))*comb(k-1, jj-1)
    # relative error
    err = abs(s - A_d[k])/abs(A_d[k])
    worst = max(worst, err)
print("mixing identity worst rel err (k<=120):", mp.nstr(worst,4))

# unambiguous phase fits: b ~ C R^-k k^-3/2 cos(k th + phi)
def fit_phase(b, k0=80, k1=120):
    rows=[]; y=[]
    for k in range(k0,k1+1):
        c = cos(k*th); s = sin(k*th)
        row=[c,s]
        for jjx in range(1,5):
            row += [c/k**jjx, s/k**jjx]
        rows.append(row); y.append(b[k]*R**k*k**mpf('1.5'))
    n=len(rows[0]); A=matrix(rows); Y=matrix(y)
    ATA=matrix(n,n); ATY=matrix(n,1)
    for i in range(n):
        for jj2 in range(n): ATA[i,jj2]=sum(A[r,i]*A[r,jj2] for r in range(len(rows)))
        ATY[i]=sum(A[r,i]*Y[r] for r in range(len(rows)))
    x=lu_solve(ATA,ATY)
    import mpmath as m
    C = sqrt(x[0]**2+x[1]**2); phi = m.atan2(x[1], x[0])
    return C, phi
bt = [A_tau[k]/factorial(k) for k in range(121)]
bd = [A_d[k]/factorial(k) for k in range(121)]
Ct, pt_ = fit_phase(bt)
Cd, pd_ = fit_phase(bd)
print("tau: C =", mp.nstr(Ct,14), " phi =", mp.nstr(pt_,12), " (theory C=",mp.nstr(sqrt(3*pi*R/8),14)," phi= pi-th/2 =", mp.nstr(pi-th/2,12),")")
print("d  : C =", mp.nstr(Cd,14), " phi =", mp.nstr(pd_,12), " (theory C=",mp.nstr(sqrt(3*pi*R/4),14)," phi= 3pi/2-th/2 =", mp.nstr(3*pi/2-th/2,12),")")
print("phi_d - phi_tau =", mp.nstr(pd_-pt_,12), " vs +pi/2 =", mp.nstr(pi/2,12), " (= -3pi/2 mod 2pi:", mp.nstr(-3*pi/2,12),")")
