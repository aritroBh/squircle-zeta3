import pickle, mpmath as mp, subprocess, time
down=pickle.load(open('gamma_down.pkl','rb')); small=pickle.load(open('gamma_small.pkl','rb'))
extra=pickle.load(open('gamma_extra.pkl','rb')); batch=pickle.load(open('gamma_batch.pkl','rb'))
def get(a1want, pools, tol=1e-6):
    best=None
    for pool in pools:
        for r in pool:
            if abs(float(r[0])-a1want)<tol: return r
    return None
pools=[small,down,extra,batch]
want=[0.00005,0.0002,0.0005,0.01,0.123,0.2,1.004406779661017,2.0]
pts=[]
for w in want:
    r=get(w,pools)
    if r: pts.append(r)
    else: print("missing", w)
print("batch points (a1):", [mp.nstr(r[0],7) for r in pts])
assert len(pts)==8
pickle.dump(pts, open('gamma_final_pts.pkl','wb'))
print("saved gamma_final_pts.pkl")
w12f = open('run_gamma_final.py').read()  # reuse script (it reads gamma_final_pts.pkl)
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_gamma_final.py'],
                      stdout=open('run_gamma_final.out','w'), stderr=subprocess.STDOUT)
open('run_gamma_final.pid','w').write(str(proc.pid))
print("launched final W12 batch pid", proc.pid)
