import time, sys, json
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps = 50
import corner_series as cs

def Mk_fft(a, Kmax, P=1024, r=mp.mpf('0.9'), verbose=False):
    """Extract M_k(a) for k=0..Kmax via Cauchy/FFT. Uses f(conj u)=conj f(u)."""
    # sample M at u_p = r e^{2pi i p/P}, p=0..P-1
    vals = [None]*P
    for p in range(P):
        th = 2*mp.pi*p/P
        u = r*mp.e**(1j*th)
        if p==0 or 2*p==P:
            vals[p] = cs.M_exact(u, a)
        elif p < P/2:
            vals[p] = cs.M_exact(u, a)
        else:
            vals[p] = mp.conj(vals[P-p])
        if verbose and p%128==0:
            print(f"  p={p}/{P} t={mp.nstr(th,4)}")
    # DFT: M_k = (1/(P r^k)) sum_p vals[p] e^{-2pi i k p/P}
    out=[]
    for k in range(Kmax+1):
        s = mp.fsum(vals[p]*mp.e**(-2j*mp.pi*k*p/P) for p in range(P))
        out.append(s/(P*r**k))
    return out

t=time.time()
a=mp.pi/2
K=40
coeffs_fft = Mk_fft(a, K, P=256, r=mp.mpf('0.9'))
print("FFT time (P=256):", round(time.time()-t,1),"s")
# compare with corner_coeffs
coeffs_ref = cs.corner_coeffs(a, 20)
for k in [0,1,2,5,10,20]:
    err = abs(coeffs_fft[k]-coeffs_ref[k])
    print(f"k={k}: fft={mp.nstr(mp.re(coeffs_fft[k]),18)} ref={mp.nstr(coeffs_ref[k],18)} |err|={mp.nstr(err,3)}")
