import os
cwd='<project>/research/08_theory/DCYA1_corner_decay_proof/'
os.chdir(cwd)
sweep = '''"""sweep_ei.py — cell sweep of E_I over [100, 10^6] with multiprocessing."""
import sys, math, json, time
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
from arb_ei import E_I_sc, iv
from multiprocessing import Pool

def make_cells():
    cells=[]
    # [100, 20000]: d(sqrt N)=0.15
    s=10.0
    while s < math.sqrt(20000):
        s2=min(s+0.15, math.sqrt(20000))
        cells.append((int(round(s*s)), int(math.ceil(s2*s2)) if s2<math.sqrt(20000) else 20000))
        s=s2
    # [20000, 10^6]: d(sqrt N)=1.0
    s=math.sqrt(20000)
    while s < 1000:
        s2=min(s+1.0, 1000)
        cells.append((int(round(s*s)), int(math.ceil(s2*s2)) if s2<1000 else 10**6))
        s=s2
    # dedupe & fix
    out=[]
    for a,b in cells:
        if b<=a: b=a+1
        if not out or a>=out[-1][1]:
            out.append((a,b))
        else:
            out[-1]=(out[-1][0],b)
    return out

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
    res=[r for r in res if r[2]==r[2]]  # drop nan
    best=max(res, key=lambda r:r[2])
    print("max E_I/sc =", best[2], "at cell", best[0], best[1])
    # print top 10
    res.sort(key=lambda r:-r[2])
    print("top 10:", [(r[0],round(r[2],2)) for r in res[:10]])
    json.dump([(r[0],r[1],r[2]) for r in res], open("ei_sweep.json","w"))
    print("time:", round(time.time()-t,1),"s")
'''
open(cwd+'sweep_ei.py','w').write(sweep)
print("sweep_ei.py written; running...")
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'sweep_ei.py'],capture_output=True,text=True,cwd=cwd,timeout=1200)
print(r.stdout)
print("STDERR:",r.stderr[-2000:] if r.stderr else "")
