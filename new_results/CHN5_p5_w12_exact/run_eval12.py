
import sys, time, pickle
sys.path.insert(0,'.')
import point_engine as PE, mzv_exact as M, kimilib
from fractions import Fraction as F
start=int(sys.argv[1]); count=int(sys.argv[2])
dd=pickle.load(open('det12_points.pkl','rb'))
sel=[tuple(x) for x in dd['sel']]
chunk=sel[start:start+count]
T=PE.table(12); d=T.nullity
def eval_point(a):
    a=[F(x) for x in a]
    g=[F(0)]*d
    for pt in (a, a[::-1]):
        vw=PE.volreg_words_point(pt,12)
        for w,pp in vw[12].items():
            rw=T.reduce_word(w)
            for b in range(d):
                if rw[b]: g[b]+=pp*rw[b]
    m0=kimilib.m0_coords(a,12)[12]
    return g, m0
out=[]
logf=open(f'eval12_{start}.log','w')
t0=time.time()
for i,a in enumerate(chunk):
    g,m0=eval_point(a)
    out.append({'pt':a,'g':g,'m0':m0})
    if (i+1)%5==0:
        logf.write(f"{i+1}/{len(chunk)} {time.time()-t0:.0f}s\n"); logf.flush()
    if (i+1)%25==0:
        pickle.dump(out, open(f'eval12_{start}.pkl','wb'))
pickle.dump(out, open(f'eval12_{start}.pkl','wb'))
logf.write(f"DONE {len(chunk)} pts {time.time()-t0:.0f}s\n"); logf.close()
print("done", start)
