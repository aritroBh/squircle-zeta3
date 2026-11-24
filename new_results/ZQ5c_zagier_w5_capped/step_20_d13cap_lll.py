import sys, json, time, pickle, itertools
from fractions import Fraction as Ff
import mpmath as mp
import flint
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
J = json.load(open('out/d13/candidates_cap.json'))
cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),13), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in J['cands']]
names=J['names']; d=len(names); N=len(cands)
Kk = json.load(open('out/d13/kernel5_cap.json'))
K=Kk['K']; comps=Kk['comps']
raw = pickle.load(open('out/d13/vals_cap_2000.pkl','rb'))
mp.mp.dps=2050
vals={k: mp.mpf(v) for k,v in raw.items()}
li5v=[vals[f'li5:{i}'] for i in range(N)]
lg=[vals[f'log:{names[i]}'] for i in range(d)]
def cv(row, comp):
    v=sum(mp.mpf(c)*li5v[i] for i,c in enumerate(row) if c)
    for ms_s, cf_s in comp.items():
        ms=eval(ms_s); cf=Ff(cf_s); pv=mp.mpf(1)
        for ii in ms: pv*=lg[ii]
        v -= mp.mpf(cf.numerator)/cf.denominator*pv
    return v
combos=[cv(r,c) for r,c in zip(K,comps)]
keepc=[j for j,v in enumerate(combos) if abs(v)>mp.mpf('1e-1900')]
print(f"d13-cap: N={N} kernel {len(K)}, EFFECTIVE combos {len(keepc)}", flush=True)
FULL = __import__('os').environ.get('FULL','0')=='1'
atoms=[]
for ms in itertools.combinations_with_replacement(range(d),3):
    v=vals['pi2']
    for i in ms: v*=lg[i]
    atoms.append((f"pi2*{ '*'.join(names[i] for i in ms) }", v))
for i in range(d): atoms.append((f"pi4*{names[i]}", vals['pi2']**2*lg[i]))
if FULL:
    atoms.append(("z3*pi2", vals['z3']*vals['pi2']))
    for a,b in itertools.combinations_with_replacement(range(d),2):
        atoms.append((f"z3*{names[a]}*{names[b]}", vals['z3']*lg[a]*lg[b]))
    atoms.append(("T3*pi2", vals['T3']*vals['pi2']))
    for a,b in itertools.combinations_with_replacement(range(d),2):
        atoms.append((f"T3*{names[a]}*{names[b]}", vals['T3']*lg[a]*lg[b]))
D5=13**5
atom_names=['T13^5','z5']+[f'K{j}' for j in keepc]+[a[0] for a in atoms]
atom_vals=[D5*vals['T'], vals['z5']]+[combos[j] for j in keepc]+[a[1] for a in atoms]
n=len(atom_vals)
tag='full' if FULL else 'd17'
print(f"d13-cap [{tag}] LLL dim {n} (combos {len(keepc)}, atoms {len(atoms)})", flush=True)
pickle.dump(dict(atom_names=atom_names, atom_vals=[mp.nstr(v,2040) for v in atom_vals], keepc=keepc),
            open(f'out/d13/atoms5_cap_{tag}.pkl','wb'))
def run(P):
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
    print(f"P={P}: LLL {dt:.0f}s, rels: {len(rels)}", flush=True)
    for mx,a,rs in rels[:5]:
        print(f"   maxcoef {mx} T-coef {a[0]} resid {mp.nstr(rs,3)} nterms {sum(1 for x in a if x)}", flush=True)
    out=[dict(maxcoef=mx, a=[str(x) for x in a], resid=mp.nstr(rs,5)) for mx,a,rs in rels[:15]]
    json.dump(dict(P=P, dt=dt, nrel=len(rels), rels=out, pool=f'cap_{tag}'), open(f'out/d13/lll_cap{tag}_P{P}.json','w'))
    return dt, [a for mx,a,rs in rels]
def norm(a):
    a=tuple(a)
    for x in a:
        if x: return tuple(-y for y in a) if x<0 else a
    return a
dt15,r15=run(1500)
dt18,r18=run(1800)
common={norm(a) for a in r15}&{norm(a) for a in r18}
print(f"1500:{len(r15)} 1800:{len(r18)} common:{len(common)}; dt {dt15:.0f}/{dt18:.0f}", flush=True)
for a in list(common)[:3]:
    print("  COMMON maxcoef", max(abs(x) for x in a), "T-coef", a[0], flush=True)
with open('progress.log','a') as f:
    f.write(f"step_20_d13cap_lll[{tag}]: dim {n}; 1500:{len(r15)} 1800:{len(r18)} common:{len(common)}\n")
