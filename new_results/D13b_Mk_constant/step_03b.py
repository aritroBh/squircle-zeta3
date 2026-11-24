
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
R = {k: d[k]*k*mp.exp(mp.sqrt(mp.pi*k))/spi for k in range(2,401)}
psi = {k: mp.sqrt(mp.pi*k) for k in range(1,402)}

def est(k):
    dk = psi[k+1]-psi[k]
    re = R[k]
    im = (R[k+1]-R[k]*mp.cos(dk))/mp.sin(dk)
    A = mp.hypot(re, im)
    th = -mp.atan2(im, re)      # theta = psi_k + phi_k
    ph = th - psi[k]
    ph = mp.fmod(ph, 2*mp.pi)
    if ph < 0: ph += 2*mp.pi
    return A, ph

print(f"{'k':>4} {'A_k est':>22} {'phi_k est':>22} {'phi_k/pi':>12}")
for k in [50,60,80,100,120,150,200,250,300,350,398]:
    A, ph = est(k)
    print(f"{k:>4} {mp.nstr(A,18):>22} {mp.nstr(ph,18):>22} {mp.nstr(ph/mp.pi,12):>12}")
