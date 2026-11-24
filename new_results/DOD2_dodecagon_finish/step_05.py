
import mpmath as mp, pickle, time
import sympy as sp
mp.mp.dps = 620
t0 = time.time()
D1 = '../DOD1_dodecagon/'
r3 = mp.sqrt(3)
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
l2,l3,lu,z2,z3 = atoms['l2'],atoms['l3'],atoms['lu'],atoms['z2'],atoms['z3']
L2x,L3x,L3o = atoms['L2'],atoms['L3'],atoms['L3o']
GD = [('3-s', 3-r3), ('s-1', r3-1), ('2s-3', 2*r3-3), ('s-2', r3-2), ('(s-1)/2', (r3-1)/2),
      ('(2s-3)/3', (2*r3-3)/3), ('2-s', 2-r3), ('(3-s)/2', (3-r3)/2), ('3-2s', 3-2*r3),
      ('4-2s', 4-2*r3), ('2-2s/3', 2-2*r3/3), ('2s-2', 2*r3-2)]
Xs = {1: 1-2*r3/3, 2: 1-r3, 5: 1-r3/3, 6: 1-r3/2}
IDX=[1,2,5,6]
extra = set()
gamS = {}
for a in IDX:
    for k in IDX:
        if a==k: continue
        g_ = sp.simplify(Xs[k]/(Xs[k]-Xs[a])); d_ = sp.simplify(g_*(1-Xs[a]))
        gamS[(a,k)] = g_
        extra.add(sp.simplify((g_-1)/g_)); extra.add(sp.simplify((d_-1)/d_))
ES = sorted(extra, key=lambda e: float(sp.N(e,10)))
En = [mp.mpf(sp.N(e,60)) for e in ES]
print("e-points:", [str(e) for e in ES])
vals=[]; names=[]
def add(nm,val): names.append(nm); vals.append(val)
add('z3',z3)
for nm,lv in [('l2',l2),('l3',l3),('lu',lu)]: add('z2'+nm,z2*lv)
for p in range(4):
    for q in range(4-p):
        for rr in range(4-p-q):
            if p+q+rr==3: add('c%d%d%d'%(p,q,rr), l2**p*l3**q*lu**rr)
add('L3x1',L3x[1]); add('L3x2',L3x[2]); add('L3o1',L3o[1]); add('L3o2',L3o[2])
for nm,lv in [('l2',l2),('l3',l3),('lu',lu)]:
    add('L2x1'+nm,L2x[1]*lv); add('L2x2'+nm,L2x[2]*lv)
for k,vv in GD: add('L3g_'+k, mp.re(mp.polylog(3, mp.mpc(vv))))
reps = ['3-2s','s-2','(2s-3)/3','2-s','(s-1)/2','2s-3']
GDd = dict(GD)
for k in reps:
    for nm,lv in [('l2',l2),('l3',l3),('lu',lu)]:
        add('L2g%s*%s'%(k,nm), mp.re(mp.polylog(2, mp.mpc(GDd[k])))*lv)
for i,e in enumerate(ES): add('L3e_%d'%i, mp.re(mp.polylog(3, mp.mpc(En[i]))))
for i,e in enumerate(ES):
    L2e = mp.re(mp.polylog(2, mp.mpc(En[i])))
    for nm,lv in [('l2',l2),('l3',l3),('lu',lu)]:
        add('L2e_%d*%s'%(i,nm), L2e*lv)
N = len(vals)
print("full basis size:", N, "t=%ds"%(time.time()-t0))
pickle.dump({'names':names,'vals':vals,'ES':[str(e) for e in ES]}, open('basis102_620.pkl','wb'))

def lll_rels(vals, C, maxc=10**6, res_thr=mp.mpf('1e-500')):
    N = len(vals)
    rows = []
    for i in range(N):
        rr_ = [0]*N + [int(mp.nint(C*vals[i]))]
        rr_[i] = 1
        rows.append(rr_)
    M = sp.Matrix(rows)
    B = M.lll(delta=sp.Rational(3,4))
    rels = []
    for i in range(B.rows):
        vc = [int(B[i,j]) for j in range(N)]
        nrm = max(abs(x) for x in vc)
        if nrm==0 or nrm > maxc: continue
        resid = abs(sum(mp.mpf(vc[q])*vals[q] for q in range(N)))
        if resid < res_thr:
            rels.append((vc, resid, nrm))
    return rels

rels = lll_rels(vals, 10**540)
print("internal relations found:", len(rels), "t=%ds"%(time.time()-t0))
# independent relations -> quotient: eliminate pivots (Gaussian elim over QQ)
import sympy as ssp
if rels:
    Rmat = ssp.Matrix([r[0] for r in rels])
    Rr, piv = Rmat.rref()
    print("relation rank:", len(piv))
    pivset = set(piv)
    keep = [i for i in range(N) if i not in pivset]
    print("independent basis size:", len(keep))
    # expression: for each relation row, pivot element = -sum_{j not pivot} Rr[row,j]/Rr[row,pivot] * elem_j
    elim = {}
    for row, pv in enumerate(piv):
        expr = {j: -Rr[row,j]/Rr[row,pv] for j in range(N) if j!=pv and Rr[row,j]!=0}
        elim[pv] = expr
    pickle.dump({'keep':keep,'elim':{k:{j:str(c) for j,c in v.items()} for k,v in elim.items()},
                 'names':names}, open('basis_quotient.pkl','wb'))
    print("saved basis_quotient.pkl")
print("elapsed", time.time()-t0)
