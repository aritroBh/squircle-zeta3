"""sweep_ei.py — cell sweep of E_I over [100, 10^6] with multiprocessing (fixed cells)."""
import sys, math, json, time
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
from arb_ei import E_I_sc, iv
from multiprocessing import Pool

def make_cells():
    cells=[]
    # region [100, 20000]: d(sqrtN)=0.15 ; [20000, 10^6]: d(sqrtN)=0.5
    for slo,shi,ds in [(10.0, math.sqrt(20000), 0.15), (math.sqrt(20000), 1000.0, 0.5)]:
        s=slo
        while s < shi-1e-12:
            s2=min(s+ds, shi)
            n1=int(math.floor(s*s)); n2=int(math.ceil(s2*s2))
            if n2<=n1: n2=n1+1
            if not cells or n1>=cells[-1][1]:
                cells.append((n1,n2))
            else:
                # overlap: shift start to previous end
                a=max(n1, cells[-1][1]); cells.append((a, max(a+1,n2)))
            s=s2
    return cells

def work(cell):
    a,b=cell
    try:
        v=float(E_I_sc(iv(a,b)))
        return (a,b,v)
    except Exception as e:
        return (a,b,float('nan'))

if __name__=="__main__":
    cells=make_cells()
    print("num cells:",len(cells))
    t=time.time()
    with Pool(8) as p:
        res=p.map(work, cells, chunksize=8)
    res=[r for r in res if r[2]==r[2]]
    best=max(res, key=lambda r:r[2])
    print("max E_I/sc =", round(best[2],3), "at cell", best[0], best[1])
    res.sort(key=lambda r:-r[2])
    print("top 12:", [(r[0],round(r[2],2)) for r in res[:12]])
    json.dump([(r[0],r[1],r[2]) for r in res], open("ei_sweep.json","w"))
    print("time:", round(time.time()-t,1),"s")
