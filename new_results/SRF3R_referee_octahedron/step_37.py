import json, time
from mpmath import mp, mpf, log, pi, sqrt, exp, factorial, sin, cos, matrix, lu_solve
t0=time.time()
mp.dps = 100
z0 = (-log(2)+pi*1j)/2
R = abs(z0); th = mp.arg(z0)
j120 = json.load(open("../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json"))
A_tau = [mpf(x) for x in j120['eta']]
b = [A_tau[k]/factorial(k) for k in range(121)]

def linfit(b, beta, k0, k1=120, ncorr=3, Rf=None, thf=None):
    Rf = Rf or R; thf = thf or th
    rows=[]; y=[]
    for k in range(k0,k1+1):
        c = cos((k-mpf('0.5'))*thf); s = sin((k-mpf('0.5'))*thf)
        row=[c, s]
        for jj in range(1,ncorr+1):
            row += [c/k**jj, s/k**jj]
        rows.append(row); y.append(b[k]*Rf**k*k**(-beta))
    n = len(rows[0])
    A = matrix(rows); Y = matrix(y)
    ATA = matrix(n,n); ATY = matrix(n,1)
    for i in range(n):
        for jj in range(n):
            ATA[i,jj] = sum(A[r,i]*A[r,jj] for r in range(len(rows)))
        ATY[i] = sum(A[r,i]*Y[r] for r in range(len(rows)))
    x = lu_solve(ATA, ATY)
    res = sqrt(sum((sum(A[r,i]*x[i] for i in range(n))-Y[r])**2 for r in range(len(rows)))/len(rows))
    return [x[i] for i in range(n)], res

# profile over beta
print("beta profile (window 50-120, 3 corrections):")
best=None
for ib in range(-1520,-1479,2):
    beta = mpf(ib)/1000
    x, res = linfit(b, beta, 50)
    flag=""
    if best is None or res<best[1]: best=(beta,res); flag=" <--"
    if ib%10==0: print(f"  beta={float(beta):.3f}: resid={mp.nstr(res,4)}{flag}")
print("best beta:", best[0], "resid", mp.nstr(best[1],4))

# finer
bf=None
for ib in range(-1510,-1489,1):
    beta = mpf(ib)/1000
    x, res = linfit(b, beta, 50)
    if bf is None or res<bf[1]: bf=(beta,res)
print("fine best beta:", bf[0], mp.nstr(bf[1],4))
# estimate with windows
for k0 in [40,50,60,70,80,90,100]:
    x, res = linfit(b, mpf('-1.5'), k0)
    print(f"k0={k0}: C1={mp.nstr(x[0],14)} C2={mp.nstr(x[1],6)} resid={mp.nstr(res,3)}  (theory C1={mp.nstr(-sqrt(3*pi/8),14)})")
