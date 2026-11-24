import pickle, numpy as np
d=pickle.load(open('signmap.pkl','rb')); grid=d['grid']; signmap=d['signmap']
# find a both-signs cell at moderate a1,a2 (near center)
cands=[(i,j) for (i,j),s in signmap.items() if len(s)>1 and 30<i<60 and 30<j<60]
print("example both-signs cells:", cands[:5])
i,j=cands[0]
print(f"chosen cell a1={grid[i]:.4g} a2={grid[j]:.4g}")
# verify with numpy all sheets (recompute carefully)
from survey_lib import roots_pos, common_pos_a5, q11_np, gcf, deg_a4
a1v,a2v=grid[i],grid[j]
cc=np.array([gcf[k](a1v,a2v) for k in range(deg_a4,-1,-1)])
print("all sheets at this (a1,a2):")
for a4r in roots_pos(cc):
    for a5r in common_pos_a5(a1v,a2v,a4r):
        if a5r>1e-7:
            q=q11_np([a1v,a2v,1.0,a4r,a5r])
            print(f"   a4={a4r:.5f} a5={a5r:.5f} Q11={q:.6g} sign={'+' if q>0 else '-'}")
