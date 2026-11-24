import numpy as np, pickle, time
from survey_lib import roots_pos, common_pos_a5, q11_np, gcf, deg_a4
# Wide survey: map sign of Q11 on V10+ over (a1,a2) in [1e-3,1e3]^2 (coarse), all sheets
t0=time.time()
N=90
grid=np.logspace(-3,3,N)
signmap={}  # (i,j) -> set of signs seen
npts=0
boundary_small_a5=[]
for i,a1v in enumerate(grid):
    for j,a2v in enumerate(grid):
        cc=np.array([gcf[k](a1v,a2v) for k in range(deg_a4,-1,-1)])
        signs=set()
        for a4r in roots_pos(cc):
            for a5r in common_pos_a5(a1v,a2v,a4r):
                if a5r>1e-6 and max(a4r,a5r)<1e8:
                    q=q11_np([a1v,a2v,1.0,a4r,a5r])
                    signs.add(1 if q>0 else (-1 if q<0 else 0)); npts+=1
        if signs: signmap[(i,j)]=signs
print(f"survey: {npts} V10+ pts, {time.time()-t0:.1f}s, {len(signmap)} (a1,a2) cells with points")
pickle.dump({'grid':grid,'signmap':signmap}, open('signmap.pkl','wb'))
# classify cells
both=sum(1 for s in signmap.values() if len(s)>1)
onlypos=sum(1 for s in signmap.values() if s=={1})
onlyneg=sum(1 for s in signmap.values() if s=={-1})
print(f"cells with BOTH signs (Q11 changes sign among sheets): {both}")
print(f"cells only positive: {onlypos}, only negative: {onlyneg}")
# where is the boundary of the V10+ projection? which cells have points
is_=[i for (i,j) in signmap]; js=[j for (i,j) in signmap]
print(f"V10+ projection covers grid index i in [{min(is_)},{max(is_)}] (a1~[{grid[min(is_)]:.3g},{grid[max(is_)]:.3g}])")
print(f"                         j in [{min(js)},{max(js)}] (a2~[{grid[min(js)]:.3g},{grid[max(js)]:.3g}])")
