
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log, nstr, cos, sin
import json
mp.dps = 80
R = {int(k): mpf(v) for k,v in json.load(open("mkmain_cache.json")).items()}

# 4-param fit: R_k = e^{-w2} N^{-1} [(C cos w2 + D sin w2) + N^{-1/2}(C2 cos w2 + D2 sin w2)], w2=sqrt(2 pi N)
def fit4(kmin,kmax):
    ks = list(range(kmin,kmax+1))
    M = mp.matrix(len(ks),4); y = mp.matrix(len(ks),1)
    for i,k in enumerate(ks):
        N = mpf(k+1); w2 = sqrt(2*pi*N); sc = exp(-w2)/N
        M[i,0] = sc*cos(w2); M[i,1] = sc*sin(w2)
        M[i,2] = sc*cos(w2)/sqrt(N); M[i,3] = sc*sin(w2)/sqrt(N)
        y[i] = R[k]
    G = M.T*M; rhs = M.T*y
    sol = mp.lu_solve(G,rhs); res = M*sol-y
    return sol, sqrt(sum(res[i]**2 for i in range(len(ks)))/len(ks))
print("4-param fit (pred: C=D=sqrt(pi/2)=1.253314, amp=sqrt(pi)=1.772454):")
for kmin,kmax in [(100,400),(150,400),(200,400),(250,400),(300,400)]:
    sol,rms = fit4(kmin,kmax)
    C,D = sol[0],sol[1]
    amp = sqrt(C*C+D*D); ph = mp.atan2(D,C)
    print("  [%3d,%3d]: C=%s D=%s amp=%s phase=%s rms=%s" % (kmin,kmax,nstr(C,9),nstr(D,9),nstr(amp,9),nstr(ph,7),nstr(rms,2)))

# phasor demodulation: Z_k = R_k e^{w2} N ; estimate local complex envelope via Z_k and Z_{k+1}
# local frequency from Z_k: use complex variable W_k = R_k e^{w2} N approx C cos+D sin -> build "analytic signal"
# simpler: fit local phase increment: solve for envelope E_k = C+iD locally from 3 consecutive k
def local_env(k):
    # solve R_j = e^{-w2_j}/N_j Re[E e^{i w2_j}] for j=k-1,k,k+1 with slowly varying E (constant over 3)
    ks = [k-1,k,k+1]
    M = mp.matrix(3,2); y = mp.matrix(3,1)
    for i,j in enumerate(ks):
        N = mpf(j+1); w2 = sqrt(2*pi*N); sc = exp(-w2)/N
        M[i,0] = sc*cos(w2); M[i,1] = sc*sin(w2); y[i] = R[j]
    G = M.T*M; rhs = M.T*y
    sol = mp.lu_solve(G,rhs)
    return sol[0]+mpc(0,1)*sol[1]
print("\nlocal envelope E_k (should -> const if frequency right):")
for k in [150, 200, 250, 300, 350]:
    E = local_env(k)
    print("  k=%3d: |E|=%s arg=%s" % (k, nstr(abs(E),9), nstr(mp.arg(E),7)))
# frequency check: detune test - refit 2-param with frequency w = sqrt(2 pi N)*(1+delta) and find min rms
def rms_detune(delta, kmin=150, kmax=400):
    ks = list(range(kmin,kmax+1))
    M = mp.matrix(len(ks),2); y = mp.matrix(len(ks),1)
    for i,k in enumerate(ks):
        N = mpf(k+1); w2 = sqrt(2*pi*N)*(1+delta); sc = exp(-sqrt(2*pi*N))/N
        M[i,0] = sc*cos(w2); M[i,1] = sc*sin(w2); y[i] = R[k]
    G = M.T*M; rhs = M.T*y
    sol = mp.lu_solve(G,rhs); res = M*sol-y
    return sqrt(sum(res[i]**2 for i in range(len(ks)))/len(ks))
print("\nfrequency detune scan (relative rms):")
base = rms_detune(mpf(0))
for d in [-0.02,-0.01,-0.005,0,0.005,0.01,0.02]:
    print("  delta=%+.3f: rms ratio=%s" % (d, nstr(rms_detune(mpf(d))/base,6)))
