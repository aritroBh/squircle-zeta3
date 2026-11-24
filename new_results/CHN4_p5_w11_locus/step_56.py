import numpy as np, pickle
from survey_lib import roots_pos, common_pos_a5, q11_np, gcf, deg_a4
# scan a2 at a1=0.01, all sheets
a1v=0.01
rows=[]
for a2v in np.linspace(0.01, 3.0, 600):
    cc=np.array([gcf[k](a1v,a2v) for k in range(deg_a4,-1,-1)])
    for a4r in roots_pos(cc):
        for a5r in common_pos_a5(a1v,a2v,a4r):
            if a5r>1e-6:
                q=q11_np([a1v,a2v,1.0,a4r,a5r])
                rows.append((a2v,a4r,a5r,q))
# cluster into sheets by (a4) continuity
print(f"total V10+ points across a2: {len(rows)}")
# For each a2, list number of points and their Q11 signs
from collections import defaultdict
by_a2=defaultdict(list)
for (a2v,a4r,a5r,q) in rows:
    by_a2[round(a2v,3)].append((a4r,q))
# print a coarse view
ks=sorted(by_a2)
for k in ks[::60]:
    info=sorted(by_a2[k])
    print(f"a2={k:.3f}: "+"; ".join(f"a4={a4:.3g},Q11={q:.3g}" for a4,q in info))
