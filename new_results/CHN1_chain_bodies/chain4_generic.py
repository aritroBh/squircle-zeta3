"""Pilot A step 5: 4-chain body with a different exponent per coordinate, u_i = 1/p_i:
V/16 = G(1+u2)G(1+u1+u3)/G(1+u1+u2+u3) * 3F2(-u4, u3, 1+u1+u3; 1+u3, 1+u1+u2+u3; 1).
Restrict to the line (u1,u2,u3,u4) = u*(a1,a2,a3,a4) and extract Taylor coefficients by FFT."""
import sys
from multiprocessing import Pool
from mpmath import mp, mpf, mpc, gamma, hyp3f2, expjpi, nstr
DPS = int(sys.argv[3]) if len(sys.argv) > 3 else 75; N = int(sys.argv[4]) if len(sys.argv) > 4 else 128; R = mpf('0.04')
A = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [3, 1, 2, 5]
def F(u):
    u1, u2, u3, u4 = [a*u for a in A]
    return gamma(1+u2)*gamma(1+u1+u3)/gamma(1+u1+u2+u3)*hyp3f2(-u4, u3, 1+u1+u3, 1+u3, 1+u1+u2+u3, 1)
def f(j):
    mp.dps = DPS
    u = R*expjpi(mpf(2*j)/N)
    v = F(u)
    return j, str(v.real), str(v.imag)
if __name__ == "__main__":
    mp.dps = 30
    # sanity: direct 2D quadrature at u = 1/7
    from mpmath import quad
    u = mpf(1)/7; u1, u2, u3, u4 = [a*u for a in A]
    d = 0 and u2*u3*quad(lambda s: quad(lambda t: s**(u2-1)*t**(u3-1)*(1-s)**u1*(1-t)**u4, [0, 1-s]), [0, 1])

    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as pool:
        res = dict((j, (a, b)) for j, a, b in pool.map(f, range(N//2 + 1)))
    mp.dps = DPS
    vals = {}
    for j in range(N):
        a, b = res[j] if j <= N//2 else res[N-j]
        v = mpc(a, b); vals[j] = v if j <= N//2 else v.conjugate()
    tag = "_".join(map(str, A))
    with open(f"chain4_generic_{tag}_d{DPS}.txt", "w") as fh:
        for k in range(13):
            ck = sum(vals[j]*expjpi(-mpf(2*j*k)/N) for j in range(N))/N/R**k
            print(k, nstr(ck.real, 50), nstr(ck.imag, 3))
            fh.write(f"{k} {nstr(ck.real, 70)}\n")
