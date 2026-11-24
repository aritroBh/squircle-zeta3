import numpy as np, pickle, time
from survey_lib import v10_points, q11_np
def survey_range(grid, a5min=1e-6, maxcoord=1e6):
    out=[]
    for a1v in grid:
        for a2v in grid:
            for a4r,a5r in v10_points(a1v,a2v):
                if a5r>a5min and max(a4r,a5r)<maxcoord:
                    q=q11_np([a1v,a2v,1.0,a4r,a5r])
                    out.append((a1v,a2v,a4r,a5r,q))
    return out
t0=time.time()
grid=np.logspace(-2,2,120)
res=survey_range(grid)
print(f"120x120 over [0.01,100]^2: {len(res)} interior V10+ pts in {time.time()-t0:.1f}s")
pickle.dump(res, open('survey_120.pkl','wb'))
qs=np.array([r[4] for r in res])
print("Q11 range:", qs.min(), "to", qs.max())
print("positive:", int(np.sum(qs>0)), "negative:", int(np.sum(qs<0)))
# well-conditioned: moderate |Q11| and moderate coords
arr=np.array(res)
mod=arr[(np.abs(arr[:,4])>1e-3)&(np.abs(arr[:,4])<1e4)]
print("\nwell-conditioned points by sign:")
poss=mod[mod[:,4]>0]; negs=mod[mod[:,4]<0]
print(f"  positive: {len(poss)}, negative: {len(negs)}")
print("  sample positive:")
for r in poss[:5]: print(f"    a1={r[0]:.4g} a2={r[1]:.4g} a4={r[2]:.4g} a5={r[3]:.4g} Q11={r[4]:.6g}")
print("  sample negative:")
for r in negs[:5]: print(f"    a1={r[0]:.4g} a2={r[1]:.4g} a4={r[2]:.4g} a5={r[3]:.4g} Q11={r[4]:.6g}")
