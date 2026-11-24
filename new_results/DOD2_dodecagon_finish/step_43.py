code = r'''
import mpmath as mp, pickle, time
import sympy as sp
mp.mp.dps = 620
t0=time.time()
D1 = '../DOD1_dodecagon/'
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
Lam = atoms['Lam']; L2x = atoms['L2']
l2,l3,lu,z2,z3 = atoms['l2'],atoms['l3'],atoms['lu'],atoms['z2'],atoms['z3']
B = pickle.load(open('basis102_620.pkl','rb'))
names, vals = B['names'], B['vals']
# known relation: Lam(1,5)-Lam(2,6)-Lam(5,1)+Lam(6,2) = (z2*(l2-l3))/2 + l2*Li2(x1) - l2*Li2(x2)
ctl = Lam[(1,5)] - Lam[(2,6)] - Lam[(5,1)] + Lam[(6,2)]
rhs = z2*(l2-l3)/2 + l2*(L2x[1]-L2x[2])
print("|ctl - rhs| =", mp.nstr(abs(ctl-rhs),4))
idx = {nm:i for i,nm in enumerate(names)}
# build the relation vector check in lattice terms:
# 2*ctl - z2l2 + z2l3 - 2*L2x1l2 + 2*L2x2l2 = 0?
relcheck = 2*ctl - vals[idx['z2l2']] + vals[idx['z2l3']] - 2*vals[idx['L2x1l2']] + 2*vals[idx['L2x2l2']]
print("|2ctl - z2l2 + z2l3 - 2L2x1l2 + 2L2x2l2| =", mp.nstr(abs(relcheck),4))
ALL = vals + [ctl]; AN = names + ['CTL']
N = len(ALL)
C = 10**540
rows = []
for i in range(N):
    rr_ = [0]*N + [int(mp.nint(C*ALL[i]))]; rr_[i] = 1
    rows.append(rr_)
print("LLL start t=%ds"%(time.time()-t0))
BM = sp.Matrix(rows).lll(delta=sp.Rational(3,4))
print("LLL done t=%ds"%(time.time()-t0))
found = []
for i in range(BM.rows):
    vc = [int(BM[i,j]) for j in range(N)]
    nrm = max(abs(x) for x in vc)
    if nrm==0 or nrm>10**6: continue
    resid = abs(sum(mp.mpf(vc[q])*ALL[q] for q in range(N)))
    if resid < mp.mpf('1e-500') and vc[N-1]!=0:
        sup = [(AN[q],vc[q]) for q in range(N) if vc[q]!=0]
        found.append((nrm,resid,sup))
print("CTL relations found:", len(found))
for nrm,resid,sup in found[:5]:
    print(" maxc=%d resid=%s"%(nrm,mp.nstr(resid,3)), sup)
print("elapsed", time.time()-t0)
'''
open('step_10.py','w').write(code)
exec(code)
