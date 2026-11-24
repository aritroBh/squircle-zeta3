import mpmath, pickle, json
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
# fit c2..c9 with the 8 points j=1..8
sel = list(range(1,9))
c = extract([taus[i] for i in sel],[Ps[i] for i in sel],P0,8)
print(f"{'k':>3} {'fitted':>28} {'predicted':>28} {'agree digits':>12}")
for j in range(8):
    k = j+2
    fit, pr = c[j], mpf(preds[k])
    d = abs(fit-pr)
    digs = -mpmath.log10(d/abs(pr)) if d>0 else 99
    print(f"{k:>3} {nstr(fit,24):>28} {nstr(pr,24):>28} {nstr(digs,5):>12}")
with open('progress.log','a') as fh:
    fh.write('step_29: squircle diamond-end c2..c9 verified: layer predictions vs perimeter fit\n')
