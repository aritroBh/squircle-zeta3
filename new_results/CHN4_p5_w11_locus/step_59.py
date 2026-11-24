w12 = r'''
import sys, time, pickle
import mpmath as mp
HERE='.'
sys.path.insert(0,HERE)
import point_engine as PE, mzv_exact as M
mp.mp.dps=150
# mpf coercion
PE.F = lambda v=0, d=None: mp.mpf(v) if d is None else mp.mpf(v)/d
def coords_mpf(a, W):
    vw=PE.volreg_words_point(a,W); vwr=PE.volreg_words_point(a[::-1],W)
    T=PE.table(W); d=T.nullity; acc=[mp.mpf(0)]*d
    for src in (vw[W],vwr[W]):
        for w,pp in src.items():
            rw=T.reduce_word(w)
            for b in range(d):
                if rw[b]: acc[b]+=pp*mp.mpf(rw[b].numerator)/rw[b].denominator
    return acc
import sympy as sp
g=pickle.load(open('gamma_pt1.pkl','rb'))
a=[g['a1'],g['a2'],mp.mpf(1),g['a4'],g['a5']]
a=[mp.mpf(x) for x in a]
print("point:", [mp.nstr(x,12) for x in a])
out={}
for W in (8,10,11,12):
    t0=time.time()
    c=coords_mpf(a,W)
    out[W]=c
    scale=max(abs(x) for x in c)
    print(f"W={W} ({time.time()-t0:.0f}s) max|coord|={mp.nstr(scale,4)}")
    for b in range(len(c)):
        mon=M.BASIS[W][b]
        lbl='*'.join(','.join(map(str,w)) for w in mon)
        print(f"   [{b}] {lbl:16s} {mp.nstr(c[b],10)}")
pickle.dump(out, open('gamma_pt1_w12.pkl','wb'))
print("saved gamma_pt1_w12.pkl")
'''
open('run_w12_gamma.py','w').write(w12)
import subprocess, time
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_w12_gamma.py'],
                      stdout=open('run_w12_gamma.out','w'), stderr=subprocess.STDOUT)
open('run_w12_gamma.pid','w').write(str(proc.pid))
print("launched W12 eval pid", proc.pid)
