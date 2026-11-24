import sys; sys.path.insert(0,'.')
import numpy as np, pickle, time
from fractions import Fraction as F
import modengine as ME
import kimilib
p = 2147483647  # 2^31-1, prime; p^2 < 2^63 safe
assert all((p-1)*(p-1) < 2**63 for _ in [0])
# K test points
pts = [(1,1,1,1,1),(2,1,3,1,2),(1,2,3,4,5),(3,1,1,2,1),(1,1,2,1,3),(2,3,1,4,1),(1,3,1,1,2)]
K=len(pts)
a_arrays = tuple(np.array([a[i]%p for a in pts],dtype=np.int64) for i in range(5))
forms = pickle.load(open('forms10.pkl','rb'))
t0=time.time()
res = ME.full_coords_mod(a_arrays, 10, p)
print(f"full_coords_mod K={K} W<=10: {time.time()-t0:.1f}s")
# compare with exact forms
bad=0
for ki,a in enumerate(pts):
    ref = kimilib.full_coords_at(a, 10, forms)
    for W in range(2,11):
        for b in range(len(ref[W])):
            exact_mod = (ref[W][b].numerator%p)*pow(ref[W][b].denominator%p,p-2,p)%p if ref[W][b] else 0
            got = int(res[W][b][ki])%p
            if exact_mod!=got: bad+=1
print("mismatches:", bad, " VALIDATION:", "PASS" if bad==0 else "FAIL")
