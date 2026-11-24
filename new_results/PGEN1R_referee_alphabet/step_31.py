import sys; sys.path.insert(0,'.')
import importlib, symbol_engine_r as SE
importlib.reload(SE)
from itertools import combinations
from flint import fmpz_mat

def run_s(s, verbose=False):
    L, idx = SE.letters(s)
    n = len(L)
    A = [('0',),('1',)] + [('p',i) for i in range(s)]      # word letters
    Pts = A + [('inf',)]                                    # cross-ratio points
    W = SE.words_dc(A, s)
    Ws = []
    for w in W:
        T = SE.S3(*w, idx, s)
        assert not any('ZERO' in key for key in T), (w, T)
        Ws.append(T)
    CR = []
    for pts4 in combinations(Pts, 4):
        CR += SE.cross_ratios(list(pts4), idx, s)
    CRset = {}
    for cr, oc in CR:
        CRset[tuple(sorted(cr.items()))] = (cr, oc)
    CRs = [SE.li3_symbol(cr, oc) for cr, oc in CRset.values()]
    Ly = SE.lyndon3(n); B = [SE.bracket_tensor(w) for w in Ly]
    # sanity: dim Q3 via pairing of all basis tensors
    Mfull = fmpz_mat([[SE.pairing_row({(0,i1,i2):0},B)][0:0] for i1 in range(0)] )
    M_cr = fmpz_mat([SE.pairing_row(T, B) for T in CRs])
    M_w  = fmpz_mat([SE.pairing_row(T, B) for T in Ws])
    M_all = fmpz_mat([[M_cr[i,j] for j in range(M_cr.ncols())] for i in range(M_cr.nrows())] +
                     [[M_w[i,j] for j in range(M_w.ncols())] for i in range(M_w.nrows())])
    # shuffle gens killed?
    D = SE.shuffle_D(n)
    bad = sum(1 for g in D if any(SE.pairing_row(g,B)))
    res = dict(s=s, n=n, nCR=len(CRset), rankCR=M_cr.rank(), nwords=len(W), rankW=M_w.rank(),
               rankUnion=M_all.rank(), dimQ3=(n**3-n)//3, Dnotkilled=bad)
    print(res, ' inclusion:', res['rankUnion']==res['rankCR'])
    return res, W, Ws, CRset, B, n, L

r1 = run_s(1)
open('progress.log','a').write('step: task A s=1 rerun after fixes\n')
