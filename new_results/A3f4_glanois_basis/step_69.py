import os, pickle
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
with open(WORK+'/pivots.pickle','rb') as f: P = pickle.load(f)
freewords = [(7,b,c) for b in (4,5,6,7) for c in range(8)]
rows_int = P['rows_int']; pivot = P['pivot_rows']; nonp = P['nonp']

# exact solve: r = c . Pivot (27x27), Fractions
Pmat = [[Fr(x) for x in rows_int[i]] for i in pivot]
def solve27(A, b):
    # gaussian elimination with partial structure, A 27x27, b 27
    n = 27
    M = [row[:] + [b[i]] for i,row in enumerate(A)]
    for col in range(n):
        # find pivot
        piv = None
        for r in range(col, n):
            if M[r][col] != 0:
                piv = r; break
        assert piv is not None, 'singular'
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
    # relation: row_k - sum c_i row_{pivot_i} = 0, as vector over the 32 words
    rel = [Fr(0)]*32
    rel[k] = Fr(1)
    for ci, pi in zip(c, pivot):
        rel[pi] = -ci
    # clear denominators
    den = reduce(lambda a,b: a*b.denominator//gcd(a,b.denominator), rel, 1)
    ints = [int(q*den) for q in rel]
    g = reduce(gcd, ints)
    ints = [a//g for a in ints]
    # make leading coeff (at k) positive
    if ints[k] < 0: ints = [-a for a in ints]
    rels.append(ints)
    print(freewords[k], 'relation: max|c| =', max(abs(a) for a in ints), 'denom cleared =', den//g)
    nonzero = [(freewords[i], a) for i,a in enumerate(ints) if a]
    print('  support:', len(nonzero), 'terms:', nonzero)
with open(WORK+'/five_relations.pickle','wb') as f:
    pickle.dump({'freewords': freewords, 'relations': rels, 'nonp': nonp}, f)
with open(WORK+'/progress.log','a') as f:
    f.write('step_26: 5 exact rational relations computed via 27x27 solves; supports printed\n')
