
import mpmath as mp
mp.mp.dps = 120

def load(fn):
    out = {}
    for line in open(fn):
        parts = line.split("\t")
        k = int(parts[0]); mid = parts[1].strip()
        out[k] = mp.mpf(mid.split(" ")[0].strip("[]"))
    return out

d = load("mk_deg9.tsv")
spi = mp.sqrt(mp.pi)

# R_k = M_k * k * exp(sqrt(pi k)) / sqrt(pi)
R = {}
for k in range(2, 401):
    R[k] = d[k]*k*mp.exp(mp.sqrt(mp.pi*k))/spi

psi = {k: mp.sqrt(mp.pi*k) for k in range(1,402)}

# analytic-signal estimator from consecutive k
def est(k):
    dk = psi[k+1]-psi[k]
    re = R[k]
    im = (R[k+1]-R[k]*mp.cos(dk))/mp.sin(dk)
    A = mp.hypot(re, im)
    ph = mp.atan2(im, re) - psi[k]
    # normalize phase to (0, 2pi]
    ph = mp.fmod(ph, 2*mp.pi)
    if ph < 0: ph += 2*mp.pi
    return A, ph

print(f"{'k':>4} {'A_k est':>22} {'phi_k est':>22}")
for k in [50,60,80,100,120,150,200,250,300,350,398]:
    A, ph = est(k)
    print(f"{k:>4} {mp.nstr(A,18):>22} {mp.nstr(ph,18):>22}")
print()
print("for reference: pi =", mp.nstr(mp.pi,20), " 3pi/4 =", mp.nstr(3*mp.pi/4,20), " 5pi/4 =", mp.nstr(5*mp.pi/4,20))
