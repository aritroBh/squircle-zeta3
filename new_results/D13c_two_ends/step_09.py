
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log, nstr, gamma, quad, cos, sin
mp.dps = 80

def load_mk(fn):
    out = {}
    for line in open(fn):
        p = line.split("\t")
        if len(p) < 2: continue
        out[int(p[0])] = mpf(p[1].strip().split(" ")[0].strip("[]"))
    return out
Mk = load_mk("../D13b_Mk_constant/mk_deg9.tsv")
def H(s): return sqrt(pi)/4 * gamma(s/2)/gamma(s/2+mpf('1.5'))
rot = exp(mpc(0,pi/4))
def Mk_main(k, dps=60):
    old = mp.dps; mp.dps = dps
    N = mpf(k+1)
    def integrand(x):
        s = rot*x
        return rot * exp(mpc(0,1)*pi*s/2) * (s/(1+s))**N * H(s)/s
    try: v = quad(integrand, [0, mpf('0.5'), 2, 10, 50, mp.inf])
    finally: mp.dps = old
    return (2/pi)*v.imag

# residual R_k = M_k - M_k^main for k=100..400
import json, os
cache = "mkmain_cache.json"
if os.path.exists(cache):
    R = {int(k): mpf(v) for k,v in json.load(open(cache)).items()}
else:
    R = {}
    for k in range(100, 401):
        R[k] = Mk[k] - Mk_main(k)
    json.dump({str(k): str(R[k]) for k in R}, open(cache,"w"))
print("residuals computed:", len(R))

# fit R_k = e^{-sqrt(2 pi N)} N^{-1} (C cos w2 + D sin w2), w2 = sqrt(2 pi N)
def fitCD(kmin, kmax):
    ks = list(range(kmin, kmax+1))
    M = mp.matrix(len(ks), 2); y = mp.matrix(len(ks),1)
    for i,k in enumerate(ks):
        N = mpf(k+1); w2 = sqrt(2*pi*N)
        sc = exp(-w2)/N
        M[i,0] = sc*cos(w2); M[i,1] = sc*sin(w2)
        y[i] = R[k]
    G = M.T*M; rhs = M.T*y
    sol = mp.lu_solve(G,rhs)
    res = M*sol-y
    rms = sqrt(sum(res[i]**2 for i in range(len(ks)))/len(ks))
    return sol[0], sol[1], rms
print("\nfit C, D over windows (prediction C=D=sqrt(pi/2)=%s):" % nstr(sqrt(pi/2),10))
for kmin,kmax in [(100,400),(150,400),(200,400),(250,400),(300,400),(200,350)]:
    C,D,rms = fitCD(kmin,kmax)
    print("  [%3d,%3d]: C=%s D=%s rms=%s" % (kmin,kmax,nstr(C,10),nstr(D,10),nstr(rms,3)))

# hypothesis test on the rate: fit R_k = e^{-a sqrtN} N^{-1}(C cos(a sqrtN)+D sin(a sqrtN)) for candidate a
print("\nrate hypothesis test (rms of fit, smaller=better):")
import math
for name, a in [("sqrt(2 pi)", sqrt(2*pi)), ("sqrt(3 pi)", sqrt(3*pi)), ("2 sqrt(pi)", 2*sqrt(pi)), ("sqrt(pi)", sqrt(pi))]:
    ks = list(range(150,401))
    M = mp.matrix(len(ks), 2); y = mp.matrix(len(ks),1)
    for i,k in enumerate(ks):
        N = mpf(k+1); w = a*sqrt(N)
        sc = exp(-w)/N
        M[i,0] = sc*cos(w); M[i,1] = sc*sin(w); y[i] = R[k]
    G = M.T*M; rhs = M.T*y
    sol = mp.lu_solve(G,rhs); res = M*sol-y
    rms = sqrt(sum(res[i]**2 for i in range(len(ks)))/len(ks))
    amp = sqrt(abs(sol[0])**2+abs(sol[1])**2)
    print("  a=%s: amp=%s rms=%s" % (name, nstr(amp,8), nstr(rms,3)))
