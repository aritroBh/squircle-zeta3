# M1 session: the single-factor claim against the ACTUAL squircle kernels (WTS kern_k.pkl, poles v in {0,a,-a,-b} -> w in {0,1,-1,-r}),
# same normalisation on both sides ([k^{d-1}] of the residue), d = 2..5, c = 1/2, 3/5.
import pickle, sys
from fractions import Fraction as F
from math import factorial
R = '<project>/research/08_theory'
sys.path.insert(0, f'{R}/WTS_weight_k_structure/theory'); import fastred  # noqa
def topc(seq, d):
    cols = [seq]
    for _ in range(d - 1): cols.append([cols[-1][i+1] - cols[-1][i] for i in range(len(cols[-1]) - 1)])
    return cols[d - 1][0]/factorial(d - 1)
SQ = {k: pickle.load(open(f'{R}/WTS_weight_k_structure/theory/kern_{k}.pkl', 'rb')) for k in range(2, 15)}
def sq_res(k, key):
    P = SQ[k][key]; return {(F(p[0]), F(p[1])): (F(c[0]), F(c[1])) for p, c in P['real'].items() if c[0] or c[1]}
VMAP = {(F(0), F(0)): 'z', (F(1), F(1)): 'p1', (F(-1), F(-1)): 'm1', (F(1), F(-1)): 'mr'}
bad = 0
for tag in ['c12', 'c35']:
    K = pickle.load(open(f'{R}/ALPHA1_angle_family_across_k/kern_{tag}.pkl', 'rb')); c = F(*K['c']); r = F(*K['r']); ks = K['kernels']
    WMAP = {F(0): 'z', F(1): 'p1', F(-1): 'm1', -r: 'mr'}
    for d in range(2, 6):
        kk = list(range(d + 1, 15)); ratios = set(); extra = []
        for j in range(d//2 + 1):
            key = (d - j, j)
            cs = [{WMAP.get(F(*p), str(p)): F(*cv[0]) for p, cv in ks[k][key]['real'].items() if F(*cv[0])} for k in kk]
            ss = [{VMAP.get(p, str(p)): v[0] for p, v in sq_res(k, key).items()} for k in kk]
            for lab in set().union(*cs, *ss):
                tc = topc([x.get(lab, F(0)) for x in cs], d); ts = topc([x.get(lab, F(0)) for x in ss], d)
                if ts == 0 and tc == 0: continue
                if ts == 0 or tc == 0: extra.append((key, lab, tc, ts)); continue
                ratios.add(tc/ts)
        ok = len(ratios) == 1 and not extra
        bad += not ok
        print(f'c={c} d={d}: corner/squircle top ratios {sorted(map(str, ratios))}; unmatched {extra}; single factor: {ok}', flush=True)
print('ALL PASS' if not bad else f'{bad} FAILURES')
