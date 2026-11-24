import sys, json, time, pickle, itertools
from fractions import Fraction as Ff
import mpmath as mp
import flint
sys.path.insert(0,'.')
import w5blib as W, w5pipe
from w5blib import QD
J, cands = w5pipe.load_cands('out/d29', 29)
names=J['names']; d=len(names); N=len(cands)
Kk = json.load(open('out/d29/kernel5.json'))
K=Kk['K']; comps=Kk['comps']
raw = pickle.load(open('out/d29/vals_2000.pkl','rb'))
mp.mp.dps=2050
vals={k: mp.mpf(v) for k,v in raw.items()}
li5v=[vals[f'li5:{i}'] for i in range(N)]
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
print(f"combos {len(combos)}, kept {len(keepc)}, exact-zero {len(combos)-len(keepc)}", flush=True)
# product atoms (weight 5): pi2*log^3, pi4*log, z3*pi2, z3*log^2, T3*pi2, T3*log^2
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
D5 = 29**5
atom_names=['T29^5','z5']+[f'K{j}' for j in keepc]+[a[0] for a in atoms]
atom_vals=[D5*vals['T'], vals['z5']]+[combos[j] for j in keepc]+[a[1] for a in atoms]
n=len(atom_vals)
print(f"LLL dim {n}  (combos {len(keepc)}, product atoms {len(atoms)})", flush=True)
pickle.dump(dict(atom_names=atom_names, atom_vals=[mp.nstr(v,2040) for v in atom_vals], keepc=keepc),
            open('out/d29/atoms5.pkl','wb'))
def run(P, thresh):
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
        if abs(res)<thresh:
            rels.append((max(abs(x) for x in a), a, res))
    rels.sort(key=lambda t:t[0])
    print(f"P={P}: LLL {dt:.0f}s, candidate rels: {len(rels)} (thresh {mp.nstr(thresh,2)})", flush=True)
    for mx,a,rs in rels[:4]:
        print(f"   maxcoef {mx} T-coef {a[0]} z5-coef {a[1]} resid {mp.nstr(rs,3)} nterms {sum(1 for x in a if x)}", flush=True)
    return dt, rels
# timing probe at 400 digits first
dt400,_ = run(400, mp.mpf('1e-350'))
print(f"[timing probe] P=400 took {dt400:.0f}s; extrapolate ~{dt400*(1500/400)**2:.0f}s @1500", flush=True)
json.dump(dict(dt400=dt400), open('out/d29/lll_timing.json','w'))
