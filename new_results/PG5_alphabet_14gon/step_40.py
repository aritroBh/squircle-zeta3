code = r'''
"""PG5 step_03: repeated-letter repair check at s=5 (and s=6), PGEN1R's method.
Every repeated-letter convergent weight-3 word reduces (shuffles + reversal) to
products and N(x,y) = 2Li3(c) - 2Li3(cX) + products, c=x/(x-y), cX=(x-1)/(x-y),
plus G(1,1,z;1)=Li3(1/(1-z)). We verify for ALL ordered pairs (x,y):
  c, cX, 1-c, 1-cX are in the CR set, and s(N) mod D lies in span s(CR)+D (rank).
"""
import sys, time, pickle
sys.path.insert(0,'.')
import symbol_engine_r as SE
from itertools import combinations
from flint import fmpz_mat
from step_02_s5_ranks import build_lyndon_index, prow

def run_repair(s):
    t0=time.time()
    L, idx = SE.letters(s)
    n = len(L)
    A = [('0',),('1',)] + [('p',i) for i in range(s)]
    Pts = A + [('inf',)]
    CR = []
    for pts4 in combinations(Pts, 4):
        CR += SE.cross_ratios(list(pts4), idx, s)
    CRset = {}
    for cr, oc in CR:
        CRset[tuple(sorted(cr.items()))] = (cr, oc)
    Ly, B, inv = build_lyndon_index(n); nB=len(B)
    CRsyms = [SE.li3_symbol(cr, oc) for cr, oc in CRset.values()]
    M_cr = fmpz_mat([prow(T, inv, nB) for T in CRsyms])
    rCR = int(M_cr.rank())
    def cr_key(v): return tuple(sorted(v.items()))
    # one-minus vector helper: 1 - z for z with vector vz: compute as ratio (differences)
    # we build c, cX, 1-c, 1-cX via ratio_vec on point differences directly
    def vec_of_ratio(num_pts, den_pts):
        st, v = SE.ratio_vec([SE.diff_vec(a,b,idx,s) for a,b in num_pts],
                             [SE.diff_vec(a,b,idx,s) for a,b in den_pts])
        assert st=='OK', (num_pts, den_pts, st)
        return v
    ZERO=('0',); ONE=('1',); INF=('inf',)
    npairs=0; nCRhit=0; sN_fail=[]
    rows_extra = []
    cases = []
    for i in range(s):
        x = ('p',i)
        ys = [ONE]+[('p',j) for j in range(s) if j!=i]
        for y in ys:
            # c = x/(x-y); cX = (x-1)/(x-y); 1-c = y/(y-x)... compute each as ratio of differences
            c   = vec_of_ratio([(x,ZERO)],[(x,y)])
            cX  = vec_of_ratio([(x,ONE)],[(x,y)])
            oc  = vec_of_ratio([(ZERO,y)],[(x,y)])     # 1-c = (0-y)/(x-y) = y/(y-x) up to sign
            ocX = vec_of_ratio([(ONE,y)],[(x,y)])      # 1-cX = (1-y)/(x-y)
            ok = all(cr_key(v) in CRset for v in (c,cX,oc,ocX))
            nCRhit += ok
            # s(N) mod D = 2 s(c) - 2 s(cX)
            sN = SE.tadd(SE.li3_symbol(c, oc), {k:-v for k,v in SE.li3_symbol(cX, ocX).items()})
            sN = {k:2*v for k,v in sN.items()}
            rows_extra.append(prow(sN, inv, nB))
            cases.append((i, y[0] if y[0]!='p' else ('p',y[1]), ok))
            npairs += 1
    # G(1,1,z) = Li3(1/(1-z)): argument vector
    rev_ok = 0
    for i in range(s):
        z = ('p',i)
        v = vec_of_ratio([(ZERO,z)],[(ONE,z)])        # 1/(1-z) = (0-z)/(1-z) up to sign
        omz = vec_of_ratio([(ONE,z)],[(ONE,z)]) if False else None
        # 1 - 1/(1-z) = z/(z-1): ratio (z-0)/(z-1)
        w = vec_of_ratio([(z,ZERO)],[(z,ONE)])
        rev_ok += (cr_key(v) in CRset and cr_key(w) in CRset)
        rows_extra.append(prow(SE.li3_symbol(v,w), inv, nB))
    M_all = fmpz_mat([M_cr[i,j] for i in range(M_cr.nrows()) for j in range(0)] and [] or
                     [[M_cr[i,j] for j in range(M_cr.ncols())] for i in range(M_cr.nrows())] + rows_extra)
    rAll = int(M_all.rank())
    res = dict(s=s, npairs=npairs, pairs_CR_ok=nCRhit, rev_cases=s, rev_ok=rev_ok,
               rankCR=rCR, rankWithN=rAll, all_in_span=(rAll==rCR), secs=round(time.time()-t0,1))
    print(res, flush=True)
    pickle.dump(res, open('repair_s%d.pkl'%s,'wb'))
    open('progress.log','a').write('PG5 step_03: repeated-letter repair s=%d: %s\n'%(s,res))
    return res

if __name__=='__main__':
    for s in [int(c) for c in (sys.argv[1] if len(sys.argv)>1 else '5')]:
        run_repair(s)
'''
open('step_03_repair.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_03_repair.py','56'],capture_output=True,text=True,timeout=1200)
print(r.stdout)
print(r.stderr[-2000:] if r.stderr else '')
