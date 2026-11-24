import sys, os, json, pickle
from fractions import Fraction as Fr
from math import comb
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
base = '<project>/research/08_theory'
tab = json.load(open(f'{base}/SUB1_subleading_table/subtop_table.json'))
K = {k: pickle.load(open(f'{base}/WTS_weight_k_structure/theory/kern_{k}.pkl','rb')) for k in range(4, 10)}
# residue(-a) at level k for key (i,j) degree d = sum_M beta_{d,M} binom(k-d-1, M-1) for k>=d+1
for (key, k) in [((2,0),4),((2,0),5),((3,0),6),((4,0),7)]:
    i,j = key; d = i+j
    betas = [Fr(x) for x in tab[str(d)][f'{key}|-a|0']['betas']]  # beta_{d,1..d}
    pred = sum(betas[M-1]*comb(k-d-1, M-2+1-1+1) for M in range(1, d+1))  # careful: beta_{d,M} with binom(k-d-1, M-1)
    pred = sum(betas[M-1]*comb(k-d-1, M-1) for M in range(1, d+1))
    rec = K[k][key]['real'].get((Fr(-1),Fr(-1)), (Fr(0),)*4)
    print(key, 'k',k, 'pred', pred, 'recorded', rec[:2], 'match:', pred==rec[0] and rec[1]==0)
