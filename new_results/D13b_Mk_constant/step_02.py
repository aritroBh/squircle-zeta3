
import mpmath as mp
mp.mp.dps = 120

def load(fn):
    out = {}
    for line in open(fn):
        parts = line.split("\t")
        k = int(parts[0])
        mid = parts[1].strip()
        # strip ball enclosure
        v = mp.mpf(mid.split(" ")[0].strip("[]"))
        out[k] = v
    return out

d8 = load("mk_deg8.tsv")
d9 = load("mk_deg9.tsv")

def sig_agree(a, b):
    if a == b: return 100
    if a == 0 or b == 0: return 0
    rel = abs(a-b)/max(abs(a),abs(b))
    if rel == 0: return 100
    return -mp.log10(rel)

print(f"{'k':>4} {'M_k (deg 9)':>26} {'sig digits deg8=deg9':>20}")
worst = 100
for k in range(50, 401, 50):
    a, b = d8[k], d9[k]
    sa = sig_agree(a,b)
    worst = min(worst, sa)
    print(f"{k:>4} {mp.nstr(b,22):>26} {mp.nstr(sa,5):>20}")
print()
# also worst over all k in [50,400]
w2 = min(sig_agree(d8[k],d9[k]) for k in range(50,401))
print("worst agreement over all k in [50,400]:", mp.nstr(w2,5), "significant digits")
w3 = min(sig_agree(d8[k],d9[k]) for k in range(1,401))
kk = min(range(1,401), key=lambda k: sig_agree(d8[k],d9[k]))
print("worst over all k in [1,400]:", mp.nstr(w3,5), "at k =", kk)
