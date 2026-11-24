# M_k from the DEFINITION: M(eps) = -F(eps)/2, F = 2*int_0^{1/2} Phi, Phi = b - b^2/(a(1+sqrt(1+z)))
# (Lemma holo form), Taylor coefficients by a discrete Cauchy integral on |eps| = r.
import sys, mpmath as mp
r = mp.mpf(sys.argv[1]); N = int(sys.argv[2]); mp.mp.dps = int(sys.argv[3]); KMAX = int(sys.argv[4])
def M(eps):
    def Phi(t):
        x, y = mp.log(t), mp.log(1-t)
        bb = mp.exp((eps-1)*y)
        b2a = mp.exp((2*eps-2)*y - (eps-1)*x)
        z = mp.exp((2-2*eps)*(x-y))
        return bb - b2a/(1+mp.sqrt(1+z))
    return -mp.quad(Phi, [0, mp.mpf(1)/8, mp.mpf(1)/2])
vals = {}
for j in range(N//2+1):
    e = r*mp.expjpi(mp.mpf(2*j)/N)
    vals[j] = M(e)
for j in range(N//2+1, N): vals[j] = mp.conj(vals[N-j])
out = []
for k in range(KMAX+1):
    s = mp.fsum(vals[j]*mp.expjpi(-mp.mpf(2*j*k)/N) for j in range(N))
    Mk = mp.re(s)/N/r**k
    out.append(Mk)
    print(k, mp.nstr(Mk, mp.mp.dps-10), flush=True)
