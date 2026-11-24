
import sys, json, pickle, time, math, random
from fractions import Fraction as Ff
import mpmath as mp
import flint
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
ZQ='../ZQ5_zagier_weight5'
D=17
J = json.load(open(f'{ZQ}/out/d17/candidates.json'))
Kk = json.load(open(f'{ZQ}/out/d17/kernel.json'))['K']
cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),17), vx=c['vx'], v1m=c['v1m']) for c in J['cands']]
names=J['names']; d=len(names); N=len(cands)
kept=[j for j in range(len(Kk)) if j not in (15,18,39,61,62)]
dropped=[15,18,39,61,62]
# products list (order as in pipe.py weight5_products)
prods = W.weight5_products(names) if hasattr(W,'weight5_products') else None
# build product list manually (ZQ5 order: log5, pi2log3, pi4log)
import itertools
prod_list=[]
for ms in itertools.combinations_with_replacement(range(d),5):
    prod_list.append(("log5",ms,"*".join(names[i] for i in ms)))
for ms in itertools.combinations_with_replacement(range(d),3):
    prod_list.append(("pi2log3",ms,"*".join([names[i] for i in ms]+["pi2"])))
for i in range(d):
    prod_list.append(("pi4log",(i,),"*".join([names[i],"pi2","pi2"])))
M=len(prod_list)
print(f"N={N} candidates, products M={M}, kept K={len(kept)}")
# names layout: [T,z5,K..,80 products]
r1150 = pickle.load(open(f'{ZQ}/planner/d17_rerun_1150.pkl','rb'))
r830 = pickle.load(open(f'{ZQ}/planner/d17_rerun_830.pkl','rb'))
assert r1150['names']==r830['names']
nm = r1150['names']
assert nm[2:60]==[f"K{j}" for j in kept]
# expand a names-space row z to atom space [z5, li5(N), prod(M)]
def expand(z):
    out=[Ff(z[1])] + [Ff(0)]*N + [Ff(int(x)) for x in z[60:]]
    for jj in range(58):
        c=z[2+jj]
        if c:
            Kr=Kk[kept[jj]]
            for i in range(N):
                if Kr[i]: out[1+i]+=Ff(int(c))*Ff(int(Kr[i]))
    return out
# relation
a1 = r1150['rels'][0][1]
c_rel = expand(a1)
# cross-check against JSON
R = json.load(open(f'{ZQ}/planner/d17_relation_planner.json'))
li5map = {k: Ff(v) for k,v in R['li5'].items()}
def xkey(x): return f"{x.a}+{x.b}*sqrt17"
ok=True; nused=0
for i,c in enumerate(cands):
    k=xkey(c['x'])
    cf = li5map.get(k, Ff(0))
    if cf!=c_rel[1+i]: ok=False; print("MISMATCH", i, k, cf, c_rel[1+i])
    if cf: nused+=1
print("relation expansion matches JSON:", ok, " used Li5:", nused, " z5:", c_rel[0], " T:", a1[0])
prod_json = {k2:Ff(v) for k2,v in R['products'].items()}
pok = all(prod_json.get(p[2],Ff(0))==c_rel[1+N+pi] for pi,p in enumerate(prod_list))
print("products match JSON:", pok, " prod terms:", sum(1 for pi in range(M) if c_rel[1+N+pi]))
# zero lattice
Z=[]
for z in [zz[1] for zz in r1150['zeros']]+[zz[1] for zz in r830['zeros']]:
    Z.append(expand(z))
for j in dropped:
    zv=[Ff(0)]+[Ff(int(x)) for x in Kk[j]]+[Ff(0)]*M
    Z.append(zv)
# integer rank + LLL reduce
Zint=[ [int(x) for x in row] for row in Z]
print("zero rows:", len(Zint), " rank:", flint.fmpz_mat(Zint).rank())
t0=time.time()
Zred = flint.fmpz_mat(Zint).lll(gram='exact')
Zr=[[int(Zred[i,j]) for j in range(Zred.ncols())] for i in range(Zred.nrows())]
Zr=[row for row in Zr if any(row)]
print(f"LLL-reduced zero basis: {len(Zr)} rows ({time.time()-t0:.0f}s)")
json.dump(dict(N=N,M=M,names=names,kept=kept,
               prod_list=[[p[0],list(p[1]),p[2]] for p in prod_list],
               c_rel=[str(x) for x in c_rel], Z=[[str(x) for x in row] for row in Zr]),
          open('d17_min_setup.json','w'))
print("saved d17_min_setup.json")
