
"""
SRF3 step_02: what singularities does the Borel transform actually have?
Take the K=80 high-precision tau-coefficients of S/S(1) (H1, 150 digits),
form b_k = a_k/k!, and study b_k * |z0|^k as a sequence: one frequency
(edge, arg(z0)=1.788) or two (edge + vertex, arg=1.353)?
Method: normalize c_k = b_k * |z0|^k, then fit
   c_k = k^be * (Ce cos(1.787953k+pe)) + k^bv * (Cv cos(1.3525..k+pv))
nonlinearly over windows, and separately single-frequency fits with free beta.
Also plot-free FFT-style diagnostic via Goertzel at the two candidate angles.
"""
import mpmath as mp, json, os
mp.mp.dps = 80
base = "<project>"
h1 = json.load(open(os.path.join(base,"research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json")))
tau = [mp.mpf(x) for x in h1['eta']]   # a_k = [tau^k] S/S(1)
dco = [mp.mpf(x) for x in h1['d']]
K = 80
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0)          # 1.6086
th_e = mp.arg(z0)    # edge argument  +1.787953
th_v = mp.pi - th_e  # vertex argument +1.352638... (z_v=(log2+ipi)/2)
print("z0 =", mp.nstr(z0,20)); print("|z0|=",mp.nstr(R,20), " A=1/R=",mp.nstr(1/R,20))
print("arg edge:", mp.nstr(th_e,15), " arg vertex:", mp.nstr(th_v,15))

b = [tau[k]/mp.factorial(k) for k in range(K+1)]
c = [b[k]*R**k for k in range(K+1)]   # normalized: should be ~ k^beta * oscillation

# Goertzel power at angle theta over window [k0,K]
def goertzel(theta, k0):
    re = sum(c[k]*mp.cos(theta*k) for k in range(k0,K+1))
    im = sum(c[k]*mp.sin(theta*k) for k in range(k0,K+1))
    n = K+1-k0
    return 2*mp.sqrt(re*re+im*im)/n

print("\nGoertzel amplitude at the two candidate angles (window k0..80):")
for k0 in [30,40,50,60]:
    print(f"  k0={k0}: edge(1.787953)={mp.nstr(goertzel(th_e,k0),6)}  vertex(1.352638)={mp.nstr(goertzel(th_v,k0),6)}")

# scan angle 0.5..2.5 for peaks
print("\nangle scan (k0=50):")
best=[]
for i in range(200):
    th = 0.6 + i*0.01
    best.append((goertzel(th,50), th))
best.sort(reverse=True)
for amp, th in best[:8]:
    print(f"   theta={mp.nstr(th,4)} amp={mp.nstr(amp,6)}")
