script = r'''
import sys, time, pickle, random
sys.path.insert(0,'.')
import point_engine as PE, mzv_exact as M
from fractions import Fraction as F

def make_f_eval(W, idx):
    T = PE.table(W)
    def f(a):
        a=[F(x) for x in a]
        acc=F(0)
        for pt in (a, a[::-1]):
            vw = PE.volreg_words_point(pt, W)
            for w,pp in vw[W].items():
                c=T.reduce_word(w)[idx]
                if c: acc+=pp*c
        return acc
    return f
f11 = make_f_eval(11,8)
_ = f11((1,2,3,4,5))  # warm cache

d=pickle.load(open('interp_points.pkl','rb'))
sel=[tuple(x) for x in d['sel']]
# verification points: distinct from sel, generic positive ints
random.seed(2024)
selset=set(sel)
ver=[]
while len(ver)<48:
    p=tuple(random.randint(1,9) for _ in range(5))
    if p not in selset and p[::-1] not in selset:
        ver.append(p); selset.add(p)
allpts = sel + ver
labels = ['fit']*len(sel)+['verify']*len(ver)
out={}
t0=time.time()
with open('f11_evals.log','w') as lg:
    for i,(p,lb) in enumerate(zip(allpts,labels)):
        v=f11(p)
        out[p]=(lb,v)
        lg.write(f"{i+1}/{len(allpts)} {lb} {p} done {time.time()-t0:.1f}s\n"); lg.flush()
pickle.dump(out, open('f11_evals.pkl','wb'))
with open('f11_evals.log','a') as lg:
    lg.write(f"ALL DONE total {time.time()-t0:.1f}s, {len(out)} points\n")
print("done", len(out))
'''
open('run_f11_evals.py','w').write(script)
print("script written")
import subprocess, os
# launch in background with nice
proc = subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_f11_evals.py'],
                        stdout=open('run_f11_evals.out','w'), stderr=subprocess.STDOUT)
open('run_f11_evals.pid','w').write(str(proc.pid))
print("launched pid", proc.pid)
