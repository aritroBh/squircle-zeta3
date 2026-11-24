import sys; sys.path.insert(0,'.')
import point_engine as PE
import time
for W in (11,12):
    t0=time.time(); T=PE.table(W); print(f"table({W}): {time.time()-t0:.1f}s nullity={T.nullity}")
t0=time.time()
r = PE.reduced_at((1,2,3,4,5), 12)
print(f"one point W<=12: {time.time()-t0:.1f}s")
print("c11 at (1,2,3,4,5):", r[11])
print("c12 (idx 7..11):", r[12][7:])
import pickle
pickle.dump({'pt':(1,2,3,4,5),'c11':r[11],'c12':r[12]}, open('test_pt_w1112.pkl','wb'))
