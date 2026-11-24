import mpmath, json
from mpmath import mp, mpf, sqrt, matrix, lu_solve, nstr
mp.dps = 90
taus = [mpf(s) for s in json.load(open('ckpt_sq2_taus.json'))]
Ps   = [mpf(s) for s in json.load(open('ckpt_sq2_Ps.json'))]
P0 = 4*sqrt(2)
def extract(taus, Ps, P0, nterms):
    n = nterms
    Q = [(P-P0)/t**2 for P,t in zip(Ps,taus)]
    A = matrix(n,n); b = matrix(n,1)
    for i,t in enumerate(taus):
        for j in range(n): A[i,j] = t**j
        b[i] = Q[i]
    return lu_solve(A,b)
preds = {2:'2.32628806654629317805252031205646842997329',3:'1.69996517519250036777674394595165573772495',
         4:'-3.15693715422070561248491036997352611509387',5:'-5.99340657187070660107725941712303200938264',
         6:'20.9114781464325920489424555553583910786688',7:'63.544518077654247431367031506734473687014',
         8:'-300.899292485393445660523299799385402718211',9:'-1287.73316060792062005836415215762506093905'}
sel = list(range(4,12))
c = extract([taus[i] for i in sel],[Ps[i] for i in sel],P0,8)
print(f"{'k':>3} {'agree digits (small-tau window)':>34}")
for j in range(8):
    k = j+2
    d = abs(c[j]-mpf(preds[k]))
    digs = -mpmath.log10(d/abs(mpf(preds[k]))) if d>0 else 99
    print(f"{k:>3} {nstr(digs,6):>34}")
# residual check of fit
Q = [(P-P0)/t**2 for P,t in zip([Ps[i] for i in sel],[taus[i] for i in sel])]
res = max(abs(sum(c[j]*taus[i]**j for j in range(8)) - Q[i]) for i in range(8))
print("max interp resid:", nstr(res,5))
