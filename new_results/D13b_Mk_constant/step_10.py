
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

print("ratio of M_k to leading natural formula  -(sqrt(pi)/N) e^{-theta} cos(theta), N=k+1:")
for k in [50,100,200,300,400]:
    N = mp.mpf(k+1); th = mp.sqrt(mp.pi*N)
    lead = -(spi/N)*mp.exp(-th)*mp.cos(th)
    print(f"  k={k}: M_k/lead = {mp.nstr(d[k]/lead, 12)}")
print()
print("note |sin theta|/|cos theta| at those k:", [mp.nstr(abs(mp.sin(mp.sqrt(mp.pi*(k+1)))/mp.cos(mp.sqrt(mp.pi*(k+1)))),3) for k in [50,100,200,300,400]])
