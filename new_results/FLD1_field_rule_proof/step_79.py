import subprocess, time
WD='<project>/research/08_theory/FLD1_field_rule_proof'
# time the word-cache build + 4-angle eval for k<=4 at N=1100, dps=70
code = r'''
import sys, pickle, time
sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps=70
from hpl import coeffs_mp
sys.path.insert(0,'<project>/research/10_frontier/G1_zeta3_any_angle/')
import corner_series
store=pickle.load(open('hpl_constructed_K5.pkl','rb'))
t0=time.time()
N=1100
words=set()
for (k,zk),terms in store.items():
    if k!=4: continue
    for (e,w,q) in terms:
        if len(w)>0: words.add(w)
wc={w: coeffs_mp(w,N) for w in words}
print("cache build", len(words), "words to N=1100:", round(time.time()-t0,1),"s")
def zeta_of(zk):
    v=mp.mpf(1)
    for i,e in enumerate(zk):
        if e: v*=mp.zeta(i+2)**e
    return v
def Mk4(c):
    tot=mp.mpf(0)
    for (kk,zk),terms in store.items():
        if kk!=4: continue
        s=mp.mpf(0)
        for (e,w,q) in terms:
            qq=mp.mpf(q.numerator)/q.denominator
            if len(w)==0: s+=qq*c**e; continue
            ce=wc[w]; sv=mp.mpf(0); xp=c
            for m in range(1,N+1): sv+=ce[m]*xp; xp*=c
            s+=qq*c**e*sv
        tot+=zeta_of(zk)*s
    return tot
for tag,a in [('2pi/3',2*mp.pi/3),('pi/2',mp.pi/2),('3pi/4',3*mp.pi/4),('pi/3',mp.pi/3)]:
    t=time.time(); c=mp.cos(a/2); Mc=corner_series.corner_coeffs(a,4)[4]
    d=abs(Mk4(c)-Mc)
    print(f'  {tag}: |diff|={mp.nstr(d,3)} (~{mp.nstr(-mp.log10(d+mp.mpf("1e-999")),4)} digits) [{time.time()-t:.0f}s]')
print("total", round(time.time()-t0,1),"s")
'''
open(WD+'/tmp_4ang.py','w').write(code)
VENV='<project>/.venv/bin/python'
t=time.time()
r=subprocess.run(['nice',VENV,'tmp_4ang.py'],capture_output=True,text=True,cwd=WD,timeout=900)
print("wall",round(time.time()-t,1),"s"); print(r.stdout); print(r.stderr[-800:] if r.returncode else '')