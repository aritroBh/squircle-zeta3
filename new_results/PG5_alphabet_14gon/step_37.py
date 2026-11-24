code = r'''
"""PG5 step_02: rank computation s=1..5 (Task 1).
Reuses PGEN1R's symbol_engine_r (exact integer arithmetic), with an inverted-index
fast pairing. Predictions: rank s(CR) mod D = 2*C(s+3,4) (140 at s=5);
word rank = s(2s^2+3s+7)/6 (60 at s=5); inclusion YES.
"""
import sys, time, pickle
sys.path.insert(0,'.')
import symbol_engine_r as SE
from itertools import combinations
from flint import fmpz_mat

def build_lyndon_index(n):
    """Lyndon basis of Lie^3 + inverted index key->[(basis_idx,coeff)]."""
    Ly = SE.lyndon3(n)
    B = [SE.bracket_tensor(w) for w in Ly]
    inv = {}
    for bi,b in enumerate(B):
        for key,c in b.items():
            inv.setdefault(key, []).append((bi,c))
    return Ly, B, inv

def prow(T, inv, nB):
    """pairing row of sparse tensor T against all Lyndon basis elements, dense list."""
    row = [0]*nB
    for key,c in T.items():
        if key in inv:
            for bi,bc in inv[key]:
                row[bi] += c*bc
    return row

def run_s(s, dim_check=True):
    t0 = time.time()
    L, idx = SE.letters(s)
    n = len(L)
    A = [('0',),('1',)] + [('p',i) for i in range(s)]
    Pts = A + [('inf',)]
    # convergent distinct-letter words
    W = SE.words_dc(A, s)
    Ws = []
    for w in W:
        T = SE.S3(*w, idx, s)
        assert not any('ZERO' in key for key in T), (w, T)
        Ws.append(T)
    # cross-ratios
    CR = []
    for pts4 in combinations(Pts, 4):
        CR += SE.cross_ratios(list(pts4), idx, s)
    CRset = {}
    for cr, oc in CR:
        CRset[tuple(sorted(cr.items()))] = (cr, oc)
    CRs = [SE.li3_symbol(cr, oc) for cr, oc in CRset.values()]
    # Lyndon pairing
    Ly, B, inv = build_lyndon_index(n)
    nB = len(B)
    rows_cr = [prow(T, inv, nB) for T in CRs]
    rows_w  = [prow(T, inv, nB) for T in Ws]
    M_cr = fmpz_mat(rows_cr)
    M_w  = fmpz_mat(rows_w)
    M_all = fmpz_mat(rows_cr + rows_w)
    # shuffle generators killed?
    bad = 0
    for g in SE.shuffle_D(n):
        if any(prow(g, inv, nB)): bad += 1
    res = dict(s=s, n=n, nCR=len(CRset), rankCR=int(M_cr.rank()), nwords=len(W),
               rankW=int(M_w.rank()), rankUnion=int(M_all.rank()),
               incl=(int(M_all.rank())==int(M_cr.rank())), dimQ3=(n**3-n)//3, Dnotkilled=bad,
               secs=round(time.time()-t0,1))
    if dim_check:
        rows = []
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    rows.append(prow({(i,j,k):1}, inv, nB))
        res['dimQ3pair'] = int(fmpz_mat(rows).rank())
    print(res, flush=True)
    return res

if __name__ == '__main__':
    import os
    which = sys.argv[1] if len(sys.argv)>1 else '12345'
    out = {}
    for s in [int(c) for c in which]:
        out[s] = run_s(s, dim_check=True)
    pickle.dump(out, open('ranks_%s.pkl'%which,'wb'))
    print('predictions: rankCR = 2*C(s+3,4):', {s: 2*len(list(combinations(range(s+3),4))) for s in out})
    print('predicted rankW = s(2s^2+3s+7)/6:', {s: s*(2*s*s+3*s+7)//6 for s in out})
    open('progress.log','a').write('PG5 step_02: ranks s=%s computed: %s\n'%(which, {s:(r['rankCR'],r['rankW'],r['incl'],r.get('dimQ3pair'),r['secs']) for s,r in out.items()}))
'''
open('step_02_s5_ranks.py','w').write(code)
print("written")
# quick timing test for s=1,2 first
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_02_s5_ranks.py','12'],capture_output=True,text=True,timeout=1200)
print(r.stdout)
print(r.stderr[-2000:] if r.stderr else '')
