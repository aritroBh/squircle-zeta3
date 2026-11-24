
import mpmath as mp, pickle, time
import sympy as sp
mp.mp.dps = 620
t0=time.time()
D1 = '../DOD1_dodecagon/'
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
X = atoms['X']; T3a = atoms['T3a']
g = pickle.load(open('gdswaps620.pkl','rb'))
SwG, SwD = g['SwG'], g['SwD']
B = pickle.load(open('basis102_620.pkl','rb'))
names, vals = B['names'], B['vals']
NB54 = 54
def test_R(key):
    i,j,k = key
    R = T3a[key] - (SwG[key]-SwD[key])/X[i]
    ALL = vals[:NB54] + [R]; AN = names[:NB54] + ['R']
    N = len(ALL)
    C = 10**540
    rows = []
    for q in range(N):
        rr_ = [0]*N + [int(mp.nint(C*ALL[q]))]; rr_[q] = 1
        rows.append(rr_)
    BM = sp.Matrix(rows).lll(delta=sp.Rational(3,4))
    hits = []
    for r in range(BM.rows):
        vc = [int(BM[r,q]) for q in range(N)]
        if vc[N-1]==0: continue
        nrm = max(abs(x) for x in vc)
        resid = abs(sum(mp.mpf(vc[q])*ALL[q] for q in range(N)))
        hits.append((nrm, resid, vc))
    hits.sort()
    print(key, "R =", mp.nstr(R,12))
    for nrm,resid,vc in hits[:6]:
        exp_res = nrm*mp.mpf(10)**(-615)
        tag = "GENUINE" if resid < exp_res*1000 else "noise"
        print("   maxc=%.2e resid=%.2e expfloor=%.1e %s" % (nrm, resid, exp_res, tag))
        if resid < exp_res*1000:
            sup = [(AN[q],vc[q]) for q in range(N) if vc[q]!=0]
            print("      ", sup[:24])
for key in [(2,5,6),(1,5,6),(2,1,6),(1,2,5)]:
    test_R(key)
print("elapsed", time.time()-t0)
