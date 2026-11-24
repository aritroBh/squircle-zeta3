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
r3 = mp.sqrt(3)
i,j,k = 2,5,6
R = T3a[(i,j,k)] - (SwG[(i,j,k)]-SwD[(i,j,k)])/X[i]
Rp = R/(1+r3)   # = R_q for this triple
ALL = vals[:54] + [Rp]; AN = names[:54] + ['Rp']
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
for resid,nrm,vc in hits[:5]:
    tf = mp.mpf(nrm)*mp.mpf(10)**-615
    print("   maxc=%s resid=%s truefloor=%s" % (mp.nstr(mp.mpf(nrm),3), mp.nstr(resid,3), mp.nstr(tf,3)))
    if resid < tf*1000:
        print("      ", [(AN[q],vc[q]) for q in range(N) if vc[q]!=0][:26])
print("elapsed", time.time()-t0)
