
import mpmath as mp
mp.mp.dps = 100

def load(fn):
    out = {}
    for line in open(fn):
        parts = line.split("\t")
        out[int(parts[0])] = mp.mpf(parts[1].strip().split(" ")[0].strip("[]"))
    return out

d = load("mk_deg9.tsv")
spi = mp.sqrt(mp.pi)

br = 1/(2*spi) + spi/4 + mp.pi**mp.mpf('1.5')/48
bi = -1/(2*spi) + spi/4 - mp.pi**mp.mpf('1.5')/48
print("beta_r =", mp.nstr(br, 20))
print("beta_i =", mp.nstr(bi, 20))
print()

# R_k := -M_k N e^th/sqrt(pi) - cos th - (br cos th + bi sin th)/sqrt(N)   ; expect (dr cos+di sin)/N + O(N^-3/2)
print(f"{'k':>4} {'N*R_k (should be smooth O(1))':>32} {'N*R_k/sqrt(N) (should vanish)':>30}")
for k in [50, 80, 120, 160, 200, 250, 300, 350, 400]:
    N = mp.mpf(k+1); th = mp.sqrt(mp.pi*N)
    Rk = -d[k]*N*mp.exp(th)/spi - mp.cos(th) - (br*mp.cos(th)+bi*mp.sin(th))/mp.sqrt(N)
    print(f"{k:>4} {mp.nstr(N*Rk, 14):>32} {mp.nstr(mp.sqrt(N)*Rk, 8):>30}")
