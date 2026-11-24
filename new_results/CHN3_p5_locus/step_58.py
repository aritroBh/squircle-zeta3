import sys; sys.path.insert(0,'.')
import time
t0=time.time()
from point_engine import reduced_at
import kimilib, pickle
forms = pickle.load(open('forms10.pkl','rb'))
pts = [(1,1,1,1,1),(2,1,3,1,2),(1,2,3,4,5),(3,1,1,2,1),(1,1,2,1,3)]
maxdiff = 0
for a in pts:
    ref = kimilib.full_coords_at(a, 10, forms)
    got = reduced_at(a, 10)
    for W in range(2,11):
        for b in range(len(ref[W])):
            d = abs(ref[W][b]-got[W][b])
            if d>maxdiff: maxdiff=d
print(f"scalar engine vs forms10, W<=10, 5 points: max|diff| = {maxdiff}  ({time.time()-t0:.1f}s)")
print("VALIDATION:", "PASS" if maxdiff==0 else "FAIL")
