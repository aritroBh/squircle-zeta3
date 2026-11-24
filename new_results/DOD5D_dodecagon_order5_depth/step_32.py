import pickle, time, sys
sys.path.insert(0, ".")
from fractions import Fraction as F
import dod5d_exact as DX

terms = pickle.load(open("E5_terms.pkl","rb"))
# quick sanity of taugam cache sizes
for (kind,k,e) in [('tau',1,5),('gam',1,5),('tau',2,2),('tau',5,1)]:
    d = DX.taugam(kind,k,e)
    print(kind,k,e,"-> entries:",len(d))
t0 = time.time()
E5X = DX.expand_E5(terms)
print("expand_E5: %.1fs, monomials: %d" % (time.time()-t0, len(E5X)))
pickle.dump(E5X, open("E5_dod_expanded.pkl","wb"))
# weight distribution
from collections import Counter
wc = Counter(sum(k[1]) for k in E5X)
print("weight distribution:", dict(sorted(wc.items())))
nc = Counter(k[0] for k in E5X)
print("t-degrees:", dict(sorted(nc.items())))
