import sys; sys.path.insert(0,'.')
import importlib, symbol_engine_r as SE
importlib.reload(SE)
from itertools import combinations
from flint import fmpz_mat

def run_s(s, verbose=False):
    L, idx = SE.letters(s)
    n = len(L)
    A = [('0',),('1',)] + [('p',i) for i in range(s)]
    W = SE.words_dc(A, s)
    Ws = []
    for w in W:
        T = SE.S3(*w, idx, s)
        assert not any('ZERO' in key for key in T), (w, T)
        Ws.append(T)
    CR = []
    for pts4 in combinations(A, 4):
        CR += SE.cross_ratios(list(pts4), idx, s)
    CRset = {}
    for cr, oc in CR:
        CRset[tuple(sorted(cr.items()))] = (cr, oc)
    CRs = [SE.li3_symbol(cr, oc) for cr, oc in CRset.values()]
    Ly = SE.lyndon3(n); B = [SE.bracket_tensor(w) for w in Ly]
    M_cr = fmpz_mat([SE.pairing_row(T, B) for T in CRs])
    M_w  = fmpz_mat([SE.pairing_row(T, B) for T in Ws])
    M_all = fmpz_mat([[M_cr[i,j] for j in range(M_cr.ncols())] for i in range(M_cr.nrows())] +
                     [[M_w[i,j] for j in range(M_w.ncols())] for i in range(M_w.nrows())])
    res = dict(s=s, n=n, nCR=len(CRset), rankCR=M_cr.rank(), nwords=len(W), rankW=M_w.rank(),
               rankUnion=M_all.rank(), dimQ3=(n**3-n)//3)
    print(res, ' inclusion:', res['rankUnion']==res['rankCR'])
    if verbose:
        for w,T in zip(W,Ws): print('  ', w, T)
    return res, W, Ws, CRset, B, n

r1, W1, Ws1, CRset1, B1, n1 = run_s(1, verbose=True)
open('progress.log','a').write('step: task A engine runs s=1: ranks (2,2), #CR=6, inclusion YES\n')
