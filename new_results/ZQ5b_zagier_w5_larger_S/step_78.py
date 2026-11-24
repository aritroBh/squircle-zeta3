import os, json
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
step = r'''
import sys, json, time, pickle
from fractions import Fraction as Ff
import mpmath as mp
import flint
sys.path.insert(0,'.')
A = pickle.load(open('out/d13/atoms5.pkl','rb'))
atom_names=A['atom_names']
mp.mp.dps=1450
atom_vals=[mp.mpf(v) for v in A['atom_vals']]
n=len(atom_vals)
def run(P):
    mp.mp.dps=P+30
    Cc=mp.power(10,P)
    B=[]
    for i,v in enumerate(atom_vals):
        row=[0]*(n+1); row[i]=1; row[n]=int(mp.nint(v*Cc)); B.append(row)
    t1=time.time(); R=flint.fmpz_mat(B).lll(gram='exact'); dt=time.time()-t1
    mp.mp.dps=1450
    rels=[]
    for i in range(n):
        a=[int(R[i,j]) for j in range(n)]
        if not any(a) or a[0]==0: continue
        res=sum(mp.mpf(x)*v for x,v in zip(a,atom_vals))
        if abs(res)<mp.mpf('1e-1000'):
            rels.append((max(abs(x) for x in a), a, mp.nstr(res,4)))
    rels.sort(key=lambda t:t[0])
    out=dict(P=P, dt=dt, nrel=len(rels), rels=[(mx, a, rs) for mx,a,rs in rels[:6]])
    json.dump(out, open(f'out/d13/lll_P{P}.json','w'))
    print(f"P={P}: LLL {dt:.0f}s, GENUINE relations (resid<1e-1000): {len(rels)}", flush=True)
    for mx,a,rs in rels[:3]:
        print(f"   maxcoef {mx} T-coef {a[0]} resid {rs} nterms {sum(1 for x in a if x)}", flush=True)
for P in (600, 900, 1100):
    run(P)
print("ALL DONE", flush=True)
'''
open('step_11b_d13_lll.py','w').write(step)
import subprocess, time
os.system('rm -f step11b.log')
proc = subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','-u','step_11b_d13_lll.py'],
                        stdout=open('step11b.log','w'), stderr=subprocess.STDOUT, cwd='.')
print("launched", proc.pid)
with open('progress.log','a') as f:
    f.write('step_10_d13_vals: values @1400d; L(chi13,5)=0.973328977797491433269523388413; Li5 2nd-way worst 6.7e-1433; 118 dup pairs worst 4.5e-1439\n')
    f.write('step_11_d13_lll: combos 263 (52 exact-zero dropped), zero-symbol atoms 72, LLL dim 285; P=400 only spurious (best 1e-373)\n')
print("logged; P=600/900/1100 running in background")
