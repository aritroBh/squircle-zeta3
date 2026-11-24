"""Pilot A step 2: Taylor coefficients of V4(u)/16 (4-chain Steinmetz body) by Cauchy/FFT on |u| = r.
V4/16 = G(1+u)G(1+2u)/G(1+3u) * 3F2(-u,u,1+2u;1+u,1+3u;1); nearest singularity u = -1/2.
Saves coefficients to chain4_coeffs.txt"""
import sys
from multiprocessing import Pool
from mpmath import mp, mpf, mpc, gamma, hyp3f2, expjpi, nstr
DPS = 75; N = 128; R = mpf('0.2')
def f(j):
    mp.dps = DPS
    u = R * expjpi(mpf(2*j)/N)
    v = gamma(1+u)*gamma(1+2*u)/gamma(1+3*u)*hyp3f2(-u, u, 1+2*u, 1+u, 1+3*u, 1)
    return j, str(v.real), str(v.imag)
if __name__ == "__main__":
    js = list(range(0, N//2 + 1))
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as pool:
        res = dict((j, (a, b)) for j, a, b in pool.map(f, js))
    mp.dps = DPS
    vals = {}
    for j in range(N):
        if j <= N//2:
            a, b = res[j]; vals[j] = mpc(a, b)
        else:
            a, b = res[N-j]; vals[j] = mpc(a, '-' + b if not b.startswith('-') else b[1:])
    K = 14
    with open("chain4_coeffs.txt", "w") as fh:
        for k in range(K+1):
            ck = sum(vals[j]*expjpi(-mpf(2*j*k)/N) for j in range(N))/N / R**k
            print(k, nstr(ck.real, 60), nstr(ck.imag, 3))
            fh.write(f"{k} {nstr(ck.real, 70)}\n")
