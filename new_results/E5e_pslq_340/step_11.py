import json, time, math
from fractions import Fraction
from flint import fmpz_mat

C = json.load(open("consts360.json"))
rel = [7680000, 0, 1175040000, 49152000, 2574432000, -501811200, -111744000, 927640800,
       6229440, 397036800, 45175200, 287114250, 381542400, 77518080, 52715520, -1240896,
       -1751040, -10005440, -2977497, -24330240, -56401920]
labels = ["X"] + C["labels"]
strs = [C["vals"]["X"]] + [C["vals"][l] for l in C["labels"]]

N = 21
SCALE = 10**300
w = [int((Fraction(s)*SCALE).numerator // (Fraction(s)*SCALE).denominator +
         (1 if (Fraction(s)*SCALE) - (Fraction(s)*SCALE).numerator // (Fraction(s)*SCALE).denominator >= Fraction(1,2) else 0))
     for s in strs]
# simpler exact rounding:
w = []
for s in strs:
    f = Fraction(s)*SCALE
    fl = f.numerator // f.denominator
    w.append(fl + (1 if f - fl >= Fraction(1,2) else 0))

M = fmpz_mat([[1 if j==i else 0 for j in range(N)] + [w[i]] for i in range(N)])
t0 = time.time()
L = M.lll(delta=0.99)
print(f"LLL done in {time.time()-t0:.1f}s")

def row(i):
    return [int(L[i,j]) for j in range(N+1)]

def norm2(r):
    return sum(x*x for x in r)

nrm_rel = math.sqrt(sum(x*x for x in rel))
print(f"||rel|| = {nrm_rel:.4e}   (Dirichlet scale 10^(300/21) ~ {10**(300/21):.2e})")
print()
found = False
for i in range(6):
    r = row(i)
    last = r[-1]
    head = r[:N]
    print(f"row{i}: norm={math.sqrt(norm2(r)):.3e}  lastcoord={last:.3e}")
    # compare with ±rel
    if head == rel or head == [-c for c in rel]:
        print("   >>> row", i, "is EXACTLY", "" if head==rel else "-", "the PSLQ relation (first 21 coords)")
        found = True
print("\nLLL recovers PSLQ relation exactly:", found)
