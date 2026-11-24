import numpy as np, pickle
from survey_lib import roots_pos, common_pos_a5, q11_np, gcf, deg_a4
a1v=0.01
# fine scan near a2 in [2.5, 3.5], all sheets
print("fine scan a1=0.01, a2 in [2.5,3.5]:")
for a2v in np.linspace(2.5, 3.5, 41):
    cc=np.array([gcf[k](a1v,a2v) for k in range(deg_a4,-1,-1)])
    pts=[]
    for a4r in roots_pos(cc):
        for a5r in common_pos_a5(a1v,a2v,a4r):
            if a5r>1e-7:
                q=q11_np([a1v,a2v,1.0,a4r,a5r])
                pts.append((a4r,a5r,q))
    s="; ".join(f"a4={a4:.4f},a5={a5:.4f},Q11={q:.5g}" for a4,a5,q in pts)
    print(f"  a2={a2v:.3f}: [{len(pts)}] {s}")
