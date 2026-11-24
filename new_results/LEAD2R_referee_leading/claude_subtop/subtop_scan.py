# M1 session scan: binomial-basis Newton coefficients beta_{d,m}(r) of the kernel residues, residue(k) = sum_m C(k-d-1, m-1) beta_{d,m},
# for m = d (top, known table) and m = d-1 (sub-top, unexplored), d = 2..9, exact from WTS kern_k (k<=18) + TWOEND1 (19..21).
import pickle, sys
from fractions import Fraction as F
from math import comb, factorial
R = '<project>/research/08_theory'
sys.path.insert(0, f'{R}/WTS_weight_k_structure/theory'); import fastred  # noqa
def ref(P):
    if isinstance(P, dict) and 'real' in P: return {(F(p[0]), F(p[1])): (F(c[0]), F(c[1])) for p, c in P['real'].items() if c[0] or c[1]}
    return {(F(r.p0), F(r.p1)): (F(c.p0), F(c.p1)) for r, c in P[0].items() if c.p0 or c.p1}
KER = {k: pickle.load(open(f'{R}/WTS_weight_k_structure/theory/kern_{k}.pkl' if k <= 18 else f'{R}/TWOEND1_weight_poles/red/kern_{k}.pkl', 'rb')) for k in range(2, 22)}
PN = {(0, 0): '0', (1, 1): 'a', (-1, -1): '-a', (-1, 1): 'b', (1, -1): '-b'}
def newton(ys):  # forward differences at k0 -> binomial coefficients beta_m = Delta^{m-1} y(k0)
    out = []; cur = list(ys)
    while cur: out.append(cur[0]); cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
    return out
for d in range(2, 10):
    ks = list(range(d+1, 22)); delta = None
    print(f'--- d={d}')
    for key in [(d-j, j) for j in range(d//2+1)]:
        seq = {k: ref(KER[k][key]) for k in ks}
        poles = sorted(set().union(*seq.values()))
        for p in poles:
            for part in (0, 1):
                nb = newton([seq[k].get(p, (F(0), F(0)))[part] for k in ks])
                top, sub = nb[d-1], nb[d-2]
                if top or sub:
                    print(f'  key={key} pole={PN.get((int(p[0]),int(p[1])),p):>2} {"rat" if part==0 else "s2 "}: top={str(top):>10}  sub={str(sub):>12}  sub/top={str(sub/top) if top else "-"}')
