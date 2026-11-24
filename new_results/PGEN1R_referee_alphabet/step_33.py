import sys; sys.path.insert(0,'.')
import importlib, symbol_engine_r as SE
importlib.reload(SE)
from itertools import combinations
from flint import fmpz_mat
import time

def run_s(s):
    t0 = time.time()
    L, idx = SE.letters(s)
    n = len(L)
    A = [('0',),('1',)] + [('p',i) for i in range(s)]
    Pts = A + [('inf',)]
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
    M_cr = fmpz_mat([SE.pairing_row(T, B) for T in CRs])
    M_w  = fmpz_mat([SE.pairing_row(T, B) for T in Ws])
    M_all = fmpz_mat([[M_cr[i,j] for j in range(M_cr.ncols())] for i in range(M_cr.nrows())] +
                     [[M_w[i,j] for j in range(M_w.ncols())] for i in range(M_w.nrows())])
    D = SE.shuffle_D(n)
    bad = sum(1 for g in D if any(SE.pairing_row(g,B)))
    # dim Q3 via full pairing
    rows = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                rows.append(SE.pairing_row({(i,j,k):1}, B))
    dimQ3pair = fmpz_mat(rows).rank()
    res = dict(s=s, n=n, nCR=len(CRset), rankCR=M_cr.rank(), nwords=len(W), rankW=M_w.rank(),
               rankUnion=M_all.rank(), incl=(M_all.rank()==M_cr.rank()), dimQ3=(n**3-n)//3,
               dimQ3pair=dimQ3pair, Dnotkilled=bad, secs=round(time.time()-t0,1))
    print(res)
    return res

out = []
for s in [1,2,3,4]:
    out.append(run_s(s))
print()
print('PGEN1 table: (1,2,2),(2,10,7),(3,30,17),(4,70,34)')
for r in out:
    print(f"s={r['s']}: nCR={r['nCR']} rankCR={r['rankCR']} rankW={r['rankW']} incl={r['incl']} dimQ3={r['dimQ3pair']}/{r['dimQ3']} Dkilled={r['Dnotkilled']==0} ({r['secs']}s)")
open('progress.log','a').write('step: task A ranks s=1..4 computed\n')
