import os
cwd = os.getcwd()
lscript = r'''import sys, json, time, pickle, itertools
from fractions import Fraction as Ff
import mpmath as mp
import flint
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
# full pool + index map
JF = json.load(open('out/d29/candidates_S7_B5.json'))
full=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),29), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in JF['cands']]
names=JF['names']; d=len(names)
idx_h = [i for i,c in enumerate(full) if c['h']<=2]
cands=[full[i] for i in idx_h]
N=len(cands)
Kk = json.load(open('out/d29/kernel5_S7h2.json'))
K=Kk['K']; comps=Kk['comps']
raw = pickle.load(open('out/d29/vals_S7_2000.pkl','rb'))
mp.mp.dps=2050
vals={k: mp.mpf(v) for k,v in raw.items()}
li5v=[vals[f'li5:{i}'] for i in idx_h]
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
keepc=[j for j,v in enumerate(combos) if abs(v)>mp.mpf('1e-1900')]
print(f"h<=2: N={N} combos {len(combos)}, kept {len(keepc)}, exact-zero {len(combos)-len(keepc)}", flush=True)
atoms=[]
for ms in itertools.combinations_with_replacement(range(d),3):
    v=vals['pi2']
    for i in ms: v*=lg[i]
    atoms.append((f"pi2*{ '*'.join(names[i] for i in ms) }", v))
for i in range(d): atoms.append((f"pi4*{names[i]}", vals['pi2']**2*lg[i]))
D5 = 29**5
atom_names=['T29^5','z5']+[f'K{j}' for j in keepc]+[a[0] for a in atoms]
atom_vals=[D5*vals['T'], vals['z5']]+[combos[j] for j in keepc]+[a[1] for a in atoms]
n=len(atom_vals)
print(f"LLL dim {n} (combos {len(keepc)}, product atoms {len(atoms)})", flush=True)
pickle.dump(dict(atom_names=atom_names, atom_vals=[mp.nstr(v,2040) for v in atom_vals], keepc=keepc),
            open('out/d29/atoms5_S7h2.pkl','wb'))
def run(P, probe=False):
    mp.mp.dps=P+30
    Cc=mp.power(10,P)
    B=[]
    for i,v in enumerate(atom_vals):
        row=[0]*(n+1); row[i]=1; row[n]=int(mp.nint(v*Cc)); B.append(row)
    t1=time.time(); R=flint.fmpz_mat(B).lll(gram='exact'); dt=time.time()-t1
    mp.mp.dps=2050
    rels=[]
    for i in range(n):
        a=[int(R[i,j]) for j in range(n)]
        if not any(a) or a[0]==0: continue
        res=sum(mp.mpf(x)*v for x,v in zip(a,atom_vals))
        if abs(res)<mp.mpf('1e-900'):
            rels.append((max(abs(x) for x in a), a, res))
    rels.sort(key=lambda t:t[0])
    print(f"P={P}: LLL {dt:.0f}s, rels(resid<1e-900): {len(rels)}", flush=True)
    for mx,a,rs in rels[:4]:
        print(f"   maxcoef {mx} T-coef {a[0]} resid {mp.nstr(rs,3)} nterms {sum(1 for x in a if x)}", flush=True)
    if not probe:
        out=[dict(maxcoef=mx, a=[str(x) for x in a], resid=mp.nstr(rs,5)) for mx,a,rs in rels[:12]]
        json.dump(dict(P=P, dt=dt, nrel=len(rels), rels=out, pool='S7h2'), open(f'out/d29/lll_S7h2_P{P}.json','w'))
    return dt, [a for mx,a,rs in rels]
dt400,_ = run(400, probe=True)
print(f"[probe] P=400 {dt400:.0f}s -> est {dt400*(1500/400)**2:.0f}s @1500", flush=True)
def norm(a):
    a=tuple(a)
    for x in a:
        if x:
            return tuple(-y for y in a) if x<0 else a
    return a
dt15, r15 = run(1500)
print(f"[time] 1500: {dt15:.0f}s", flush=True)
dt18, r18 = run(1800)
print(f"[time] 1800: {dt18:.0f}s", flush=True)
common = {norm(a) for a in r15} & {norm(a) for a in r18}
print(f"1500: {len(r15)} rels, 1800: {len(r18)} rels, common: {len(common)}", flush=True)
for a in list(common)[:3]:
    print("  COMMON:", a[:6], "maxcoef", max(abs(x) for x in a), flush=True)
with open('progress.log','a') as f:
    f.write(f"step_16_d29S7h2_lll: dim {n}; 1500: {len(r15)} rels ({dt15:.0f}s), 1800: {len(r18)} ({dt18:.0f}s), common {len(common)}\n")
'''
open('step_16_d29S7_lll.py','w').write(lscript)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_16_d29S7_lll.py'],
                   capture_output=True, text=True, cwd=cwd, timeout=2700)
print("elapsed", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[-1500:] if r.stderr else "")
