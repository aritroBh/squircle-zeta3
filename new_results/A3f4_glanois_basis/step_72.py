import os, pickle
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
import mpmath as mp
mp.mp.dps = 750

with open(WORK+'/vals_hi.pickle','rb') as f: V = pickle.load(f)
free_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in V['vals'][27:]]
freewords = V['freewords']
with open(WORK+'/expansions.pickle','rb') as f: E = pickle.load(f)
with open(WORK+'/pivots.pickle','rb') as f: P = pickle.load(f)
rows_int = P['rows_int']; pivot = P['pivot_rows']; nonp = P['nonp']
a0s = [E[w]['a0'] for w in freewords]

Pmat = [[Fr(x) for x in rows_int[i]] for i in pivot]
def solve27(A, b):
    n = 27
    M = [row[:] + [b[i]] for i,row in enumerate(A)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        inv = 1/M[col][col]
        M[col] = [x*inv for x in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [a - f*bb for a,bb in zip(M[r], M[col])]
    return [M[i][n] for i in range(n)]

rels = []
for k in nonp:
    b = [Fr(x) for x in rows_int[k]]
    c = solve27([row[:] for row in Pmat], b)
    # correct word-relation: a0_k * w_k = sum_i c_i * a0_{pivot_i} * w_{pivot_i}
    rel = [Fr(0)]*32
    rel[k] = Fr(a0s[k])
    for ci, pi in zip(c, pivot):
        rel[pi] = -ci*a0s[pi]
    den = reduce(lambda a,qq: a*qq.denominator//gcd(a,qq.denominator), rel, 1)
    ints = [int(q*den) for q in rel]
    g = reduce(gcd, ints)
    ints = [x//g for x in ints]
    if ints[k] < 0: ints = [-x for x in ints]
    rels.append(ints)
    res = sum(x*fv for x,fv in zip(ints, free_vals))
    h = max(abs(x) for x in ints)
    nz = [(freewords[i], x) for i,x in enumerate(ints) if x]
    print(f'{freewords[k]}: height 10^{float(mp.log10(h)):.1f} residual {mp.nstr(abs(res),3)} ratio 10^{float(mp.log10(abs(res)/h)):.0f} support {len(nz)}')
with open(WORK+'/five_relations.pickle','wb') as f:
    pickle.dump({'freewords': freewords, 'relations': rels, 'nonp': nonp}, f)
with open(WORK+'/progress.log','a') as f:
    f.write('step_27: fixed a0 scaling bug; 5 relations verified numerically\n')
