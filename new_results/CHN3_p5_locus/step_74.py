src = open('modengine.py').read()
src = src.replace("        q[m]=((((-1)**m)%p)*invp[m]%p)*s%p",
                  "        q[m]={(m,): ((((-1)**m)%p)*invp[m]%p)*s%p}")
open('modengine.py','w').write(src)
import importlib, sys; sys.path.insert(0,'.')
import numpy as np, pickle, time
from fractions import Fraction as F
import modengine as ME
importlib.reload(ME)
import kimilib
p = 2147483647
pts = [(1,1,1,1,1),(2,1,3,1,2),(1,2,3,4,5),(3,1,1,2,1),(1,1,2,1,3),(2,3,1,4,1),(1,3,1,1,2)]
K=len(pts)
a_arrays = tuple(np.array([a[i]%p for a in pts],dtype=np.int64) for i in range(5))
forms = pickle.load(open('forms10.pkl','rb'))
t0=time.time()
res = ME.full_coords_mod(a_arrays, 10, p)
print(f"K={K} W<=10: {time.time()-t0:.1f}s")
bad=0
for ki,a in enumerate(pts):
    ref = kimilib.full_coords_at(a, 10, forms)
    for W in range(2,11):
        for b in range(len(ref[W])):
            exact_mod = (ref[W][b].numerator%p)*pow(ref[W][b].denominator%p,p-2,p)%p if ref[W][b] else 0
            if exact_mod != int(res[W][b][ki])%p: bad+=1
print("mismatches:", bad, " VALIDATION:", "PASS" if bad==0 else "FAIL")
