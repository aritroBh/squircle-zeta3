"""sweep fine: full range [100, 10^6], d(sqrtN)=0.05, find global max E_I/sc."""
import sys, math, json, time
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
from arb_ei import E_I_sc, iv
from multiprocessing import Pool
def make_cells():
    cells=[]; s=10.0
    while s < 1000:
        s2=min(s+0.05, 1000.0)
        n1=int(math.floor(s*s)); n2=int(math.ceil(s2*s2))
        if n2<=n1: n2=n1+1
        if not cells or n1>=cells[-1][1]: cells.append((n1,n2))
        else: cells.append((max(n1,cells[-1][1]), max(n2,cells[-1][1]+1)))
        s=s2
    return cells
def work(cell):
    a,b=cell
    try: return (a,b,float(E_I_sc(iv(a,b))))
    except Exception: return (a,b,float('nan'))
if __name__=="__main__":
    cells=make_cells(); print("num cells:",len(cells))
    t=time.time()
    with Pool(8) as p:
        res=p.map(work, cells, chunksize=16)
    res=[r for r in res if r[2]==r[2]]
    res.sort(key=lambda r:-r[2])
    print("top 15:", [(r[0],round(r[2],2)) for r in res[:15]])
    # profile summary at decades
    res.sort()
    for target in [100,400,1600,6400,25600,100000,400000,1000000]:
        cand=[r for r in res if r[0]<=target<=r[1]]
        if cand: print(f"N~{target}: {cand[0][2]:.2f}")
    json.dump([(r[0],r[1],r[2]) for r in res], open("ei_sweep_fine.json","w"))
    print("time:", round(time.time()-t,1),"s")
