import os, json
from fractions import Fraction as Ff
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
step = r'''
import sys, json, time, pickle, itertools
from fractions import Fraction as Ff
import mpmath as mp
import flint
sys.path.insert(0,'.')
import w5blib as W, w5pipe
from w5blib import QD
J = json.load(open('out/d13/candidates.json'))
allc=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),13), vx=c['vx'], v1m=c['v1m']) for c in J['cands']]
Sset=set(c['x'] for c in allc)
keep=[]
for c in allc:
    red=False
    for y in Sset:
        if y*y==c['x'] and -y in Sset: red=True; break
    if not red: keep.append(c)
cands=keep
names=J['names']; d=len(names); N=len(cands)
print(f"reduced pool N={N} (from {len(allc)})", flush=True)
K, comps, sym_dicts, rank = w5pipe.kernel_and_comps(cands, d, 'd13r')
json.dump(dict(K=K, comps=comps, names=names, d=d, rank=rank, reduced=True),
          open('out/d13/kernel5_red.json','w'))
print("saved kernel5_red.json", flush=True)
# build atoms
raw = pickle.load(open('out/d13/vals_1400.pkl','rb'))
mp.mp.dps=1450
vals={k: mp.mpf(v) for k,v in raw.items()}
# need li5 for kept candidates: recompute indices via x match
xkey=lambda x: x
oldmap={c['x']: i for i,c in enumerate([dict(x=QD(Ff(cc['x'][0]),Ff(cc['x'][1]),13)) for cc in J['cands']])}
li5v=[vals[f'li5:{oldmap[c["x"]]}'] for c in cands]
lg=[vals[f'log:{names[i]}'] for i in range(d)]
def combo_val(row, comp):
    v=sum(mp.mpf(c)*li5v[i] for i,c in enumerate(row) if c)
    for ms_s, cf_s in comp.items():
        ms=eval(ms_s); cf=Ff(cf_s)
        pv=mp.mpf(1)
        for ii in ms: pv*=lg[ii]
        v -= mp.mpf(cf.numerator)/cf.denominator*pv
    return v
combos=[combo_val(row,comp) for row,comp in zip(K,comps)]
keepc=[j for j,v in enumerate(combos) if abs(v)>mp.mpf('1e-1300')]
print(f"combos {len(combos)}, kept {len(keepc)}, exact-zero {len(combos)-len(keepc)}", flush=True)
atoms=[]
for ms in itertools.combinations_with_replacement(range(d),3):
    v=vals['pi2']
    for i in ms: v*=lg[i]
    atoms.append((f"pi2*{ '*'.join(names[i] for i in ms) }", v))
for i in range(d): atoms.append((f"pi4*{names[i]}", vals['pi2']**2*lg[i]))
atoms.append(("z3*pi2", vals['z3']*vals['pi2']))
for a,b in itertools.combinations_with_replacement(range(d),2):
    atoms.append((f"z3*{names[a]}*{names[b]}", vals['z3']*lg[a]*lg[b]))
atoms.append(("T3*pi2", vals['T3']*vals['pi2']))
for a,b in itertools.combinations_with_replacement(range(d),2):
    atoms.append((f"T3*{names[a]}*{names[b]}", vals['T3']*lg[a]*lg[b]))
atom_names=['T','z5']+[f'K{j}' for j in keepc]+[a[0] for a in atoms]
atom_vals=[vals['T'], vals['z5']]+[combos[j] for j in keepc]+[a[1] for a in atoms]
print(f"LLL dim {len(atom_vals)}", flush=True)
pickle.dump(dict(atom_names=atom_names, atom_vals=[mp.nstr(v,1440) for v in atom_vals]), open('out/d13/atoms5_red.pkl','wb'))
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
    json.dump(dict(P=P,dt=dt,nrel=len(rels),rels=[(mx,a,rs) for mx,a,rs in rels[:6]]), open(f'out/d13/lll_red_P{P}.json','w'))
    print(f"P={P}: LLL {dt:.0f}s, GENUINE rels: {len(rels)}", flush=True)
    for mx,a,rs in rels[:3]:
        print(f"   maxcoef {mx} T-coef {a[0]} resid {rs} nterms {sum(1 for x in a if x)}", flush=True)
for P in (900,1100):
    run(P)
print("ALL DONE", flush=True)
'''
open('step_12_d13_red_lll.py','w').write(step)
import subprocess, time
os.system('rm -f step12.log')
proc = subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','-u','step_12_d13_red_lll.py'],
                        stdout=open('step12.log','w'), stderr=subprocess.STDOUT, cwd='.')
print("launched", proc.pid)
time.sleep(240)
print(open('step12.log').read())
print("running:", subprocess.run(['pgrep','-f','step_12_d13_red_lll'],capture_output=True,text=True).stdout.strip() != '')
