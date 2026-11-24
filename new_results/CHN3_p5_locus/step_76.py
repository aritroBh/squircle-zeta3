import sys; sys.path.insert(0,'.')
import numpy as np, pickle
import modengine as ME
import kimilib
p=2147483647
pts = [(1,1,1,1,1),(2,1,3,1,2),(1,2,3,4,5),(3,1,1,2,1),(1,1,2,1,3),(2,3,1,4,1),(1,3,1,1,2)]
a_arrays = tuple(np.array([a[i]%p for a in pts],dtype=np.int64) for i in range(5))
forms = pickle.load(open('forms10.pkl','rb'))
res = ME.full_coords_mod(a_arrays, 10, p)
for W in range(2,11):
    badw=0; first=None
    for ki,a in enumerate(pts):
        ref = kimilib.full_coords_at(a, 10, forms)
        for b in range(len(ref[W])):
            exact_mod = (ref[W][b].numerator%p)*pow(ref[W][b].denominator%p,p-2,p)%p if ref[W][b] else 0
            if exact_mod != int(res[W][b][ki])%p:
                badw+=1
                if first is None: first=(ki,a,b,int(res[W][b][ki])%p,exact_mod)
    print(f"W={W}: {badw} mismatches", " first:", first)
