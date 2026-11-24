import sys; sys.path.insert(0,'.')
import symbol_engine_r as SE
from itertools import combinations
from flint import fmpz_mat

s = 1
L, idx = SE.letters(s)
n = len(L)
print('s=%d letters n=%d:'%(s,n), L)
A = [('0',),('1',)] + [('p',i) for i in range(s)]
W = SE.words_dc(A, s)
print('convergent distinct-letter words (#%d):'%len(W), W)
Ws = []
for w in W:
    T = SE.S3(*w, idx, s)
    assert not any('ZERO' in key for key in T), (w, T)
    Ws.append(T)
    print(w, T)
CR = []
for pts4 in combinations(A, 4):
    CR += SE.cross_ratios(list(pts4), idx, s)
CRset = {}
for cr, oc in CR:
    CRset[tuple(sorted(cr.items()))] = (cr, oc)
print('#CR distinct:', len(CRset), ' (PGEN1 table: 6)')
CRs = [SE.li3_symbol(cr, oc) for cr, oc in CRset.values()]
Ly = SE.lyndon3(n); B = [SE.bracket_tensor(w) for w in Ly]
print('dim Lie^3 =', len(Ly), ' Witt check:', (n**3-n)//3)
full = fmpz_mat([[SE.pairing_row({(i,j,k):1}, B) for k in range(n) for j in range(n)] for i in range(n)])
print('dim Q3 via pairing =', full.rank(), '(expect', (n**3-n)//3,')')
M_cr = fmpz_mat([SE.pairing_row(T, B) for T in CRs])
M_w  = fmpz_mat([SE.pairing_row(T, B) for T in Ws])
M_all = fmpz_mat([[M_cr[i,j] for j in range(M_cr.ncols())] for i in range(M_cr.nrows())] +
                 [[M_w[i,j] for j in range(M_w.ncols())] for i in range(M_w.nrows())])
print('rank s(CR) mod D =', M_cr.rank(), ' (PGEN1: 2)')
print('rank words mod D =', M_w.rank(), ' (PGEN1: 2)')
print('rank union =', M_all.rank(), ' inclusion:', M_all.rank()==M_cr.rank())
D = SE.shuffle_D(n)
bad = sum(1 for g in D if any(SE.pairing_row(g,B)))
print('D gens NOT killed by pairing:', bad, 'of', len(D))
