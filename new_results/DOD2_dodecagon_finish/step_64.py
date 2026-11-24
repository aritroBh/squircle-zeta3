code = r'''
import mpmath as mp, pickle, time
import sympy as sp
mp.mp.dps = 620
t0=time.time()
B = pickle.load(open('basis102_620.pkl','rb'))
names, vals = B['names'], B['vals']
Xv = pickle.load(open('Xi_verified_620.pkl','rb'))
targets = [('XiD_p', Xv['XiD_p']), ('XiD_q', Xv['XiD_q'])]
def lll_target(tname, tv):
    ALL = vals + [tv]; AN = names + [tname]
    N = len(ALL)
    C = 10**540
    rows = []
    for q in range(N):
        rr_ = [0]*N + [int(mp.nint(C*ALL[q]))]; rr_[q] = 1
        rows.append(rr_)
    BM = sp.Matrix(rows).lll(delta=sp.Rational(3,4))
    hits=[]
    for r in range(BM.rows):
        vc = [int(BM[r,q]) for q in range(N)]
        if vc[N-1]==0: continue
        nrm = max(abs(x) for x in vc)
        resid = abs(sum(mp.mpf(vc[q])*ALL[q] for q in range(N)))
        hits.append((resid, nrm, vc))
    hits.sort()
    print(tname, "t=%ds: %d target-involving vectors"%(time.time()-t0, len(hits)))
    genuine = []
    for resid,nrm,vc in hits[:6]:
        tf = mp.mpf(nrm)*mp.mpf(10)**-615
        genu = resid < tf*1000
        print("   maxc=%s resid=%s floor=%s %s" % (mp.nstr(mp.mpf(nrm),3), mp.nstr(resid,3), mp.nstr(tf,3), "GENUINE" if genu else "noise"))
        if genu:
            genuine.append((vc,resid,nrm))
            print("      ", [(AN[q],vc[q]) for q in range(N) if vc[q]!=0][:28])
    return genuine
for tname, tv in targets:
    lll_target(tname, tv)
print("elapsed", time.time()-t0)
'''
open('step_17.py','w').write(code)
exec(code)
