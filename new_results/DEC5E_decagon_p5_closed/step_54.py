import sys
sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 200
phi = (1+mp.sqrt(5))/2
def Li32(x, y, D=200):
    mp.mp.dps = D
    s = mp.mpf(0); k = 0
    # sum_{n>m>=1} x^n y^m / (n^3 m^2)
    N = 600
    xn = x
    for n in range(1, N+1):
        if n > 1: xn *= x
        ym = y
        for m in range(1, n):
            if m > 1: ym *= y
            s += xn*ym/(n**3 * m**2)
    return s
cands = [(-1/phi, 1), (1/phi**2, 1), (-1/phi, -1/phi), (1/phi**2, -phi)]
vals = [Li32(x,y) for (x,y) in cands]
for (x,y), v in zip(cands, vals):
    print(f"Li32({mp.nstr(x,6)}, {mp.nstr(y,6)}) = {mp.nstr(v, 20)}")
# save at 800 digits for later
vals800 = []
for (x,y) in cands:
    vals800.append(Li32(x,y,800))
import pickle
pickle.dump({'args': [(mp.nstr(x,780), mp.nstr(y,780)) for x,y in cands],
             'vals': [mp.nstr(v,780) for v in vals800]}, open('li32_gold_800.pkl','wb'))
print("saved li32_gold_800.pkl")
