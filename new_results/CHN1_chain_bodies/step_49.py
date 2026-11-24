import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
import pickle, sympy as sp
cache=pickle.load(open('locus_logV.pkl','rb'))
for m in range(2,12):
    p1=cache['Pi1']['prods'][m]; p2=cache['Pi2']['prods'][m]
    keys=set(p1)|set(p2)
    diffs=[(mon, sp.expand(p1.get(mon,0)-p2.get(mon,0))) for mon in keys if sp.expand(p1.get(mon,0)-p2.get(mon,0))!=0]
    if diffs:
        print('m=%d DIFFERS:'%m, [(str(mon),str(d)) for mon,d in diffs])
    else:
        print('m=%d identical'%m)
