import sys, time, itertools, pickle
sys.path.insert(0, ".")
from oct9d_common import *
t0=time.time()

# Candidate pools
poolA = pool([1,2,3], (2,3,5), SMALL235)          # edges 1-3
poolB = pool([0,1,2,3], (2,3,5,7), SMALL235+SMALL7)  # edge 0
print("poolA ({2,3,5}):", len(poolA), "candidates")
print(poolA)
print("poolB ({2,3,5,7}):", len(poolB), "candidates")
print(poolB)
pickle.dump({'poolA': poolA, 'poolB': poolB}, open("pools.pkl","wb"))
print(f"({time.time()-t0:.1f}s)")
