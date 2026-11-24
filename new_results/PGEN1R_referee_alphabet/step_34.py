import sys; sys.path.insert(0,'.')
import importlib, symbol_engine_r as SE
importlib.reload(SE)
from itertools import permutations
from flint import fmpz_mat
import random

# --- 1. shuffle homomorphism: S(G_u) shuffle S(G_v) == sum over u sh v of S(G_w), exact ---
def shuffle_words(u, v):
    out = []
    def rec(i, j, pre):
        if i==len(u) and j==len(v): out.append(pre); return
        if i<len(u): rec(i+1, j, pre+(u[i],))
        if j<len(v): rec(i, j+1, pre+(v[j],))
    rec(0,0,())
    return out

def tshuffle1(S1, S2p):
    """shuffle of a 1-tensor {(i,):c} with 2-tensor {(j,k):c}."""
    out = {}
    for (i,),c1 in S1.items():
        for (j,k),c2 in S2p.items():
            for key in [(i,j,k),(j,i,k),(j,k,i)]:
                out[key] = out.get(key,0)+c1*c2
                if out[key]==0: del out[key]
    return out

random.seed(1)
s = 2
L, idx = SE.letters(s)
A = [('0',),('1',)] + [('p',i) for i in range(s)]
ok = True
for trial in range(40):
    a = random.choice(A[1:])         # 1-letter word (avoid 0: divergent)
    bc = tuple(random.choice(A) for _ in range(2))
    if bc[0]==bc[1]: continue
    # word G(a) shuffle G(bc): words of length 3; symbol via S3 needs distinct letters, so only
    # test when all shuffle results have distinct letters
    ws = shuffle_words((a,), bc)
    if any(len(set(w))<3 for w in ws): continue
    lhs = tshuffle1({(k,):v for k,v in SE.s1_vec(a, idx, s).items()}, SE.S2_pair(*bc, idx, s))
    rhs = SE.tadd(*[SE.S3(*w, idx, s) for w in ws])
    if lhs != rhs:
        ok = False; print('SHUFFLE FAIL', a, bc); break
print('shuffle-homomorphism checks (exact):', 'PASS' if ok else 'FAIL')

# --- 2. formal shuffle reduction: every convergent word reduces to distinct-letter convergent
#        words modulo shuffle relations (exact linear algebra over Q) ---
def shuffle_reduction_test(s):
    L, idx = SE.letters(s)
    A = [('0',),('1',)] + [('p',i) for i in range(s)]
    W_all = list(permutations(A,2))  # placeholder
    import itertools
    W_all = list(itertools.product(A, repeat=3))
    wid = {w:k for k,w in enumerate(W_all)}
    Wdc = [w for w in W_all if len(set(w))==3 and w[0]!=('0',) and w[2]!=('1',)]
    Wconv = [w for w in W_all if w[0]!=('0',) and w[2]!=('1',)]
    rows = []
    # unit rows for distinct-letter convergent words
    for w in Wdc:
        r = [0]*len(W_all); r[wid[w]] = 1; rows.append(r)
    # shuffle relations: for u (len1), v (len2): sum over u sh v = 0  (mod products)
    for u in A:
        for v in itertools.product(A, repeat=2):
            r = [0]*len(W_all)
            for w in shuffle_words((u,), v):
                r[wid[w]] += 1
            rows.append(r)
    M = fmpz_mat(rows)
    base_rank = M.rank()
    # add each convergent word; rank must not grow
    okall = True
    extra = [w for w in Wconv if w not in set(Wdc)]
    Mextra = fmpz_mat(rows + [[1 if k==wid[w] else 0 for k in range(len(W_all))] for w in extra])
    okall = (Mextra.rank() == base_rank)
    print(f's={s}: #words^3={len(W_all)}, #dc={len(Wdc)}, #conv={len(Wconv)}, rank rows={base_rank}, '
          f'rank with all convergent={Mextra.rank()}, reduction OK: {okall}')
    return okall

for s in [1,2,3,4]:
    shuffle_reduction_test(s)
open('progress.log','a').write('step: task A shuffle homomorphism + formal shuffle reduction tests\n')
