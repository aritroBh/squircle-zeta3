
import mpmath as mp, pickle, time
import sympy as sp
mp.mp.dps = 150
t0 = time.time()
r3 = mp.sqrt(3)
X = {1: 1-2*r3/3, 2: 1-r3, 5: 1-r3/3, 6: 1-r3/2}
IDX=[1,2,5,6]
l2,l3,lu,z2,z3 = mp.log(2),mp.log(3),mp.log(2-r3),mp.pi**2/6,mp.zeta(3)
L2x = {i: mp.polylog(2,X[i]) for i in IDX}
L3x = {i: mp.polylog(3,X[i]) for i in IDX}
L3o = {1: mp.polylog(3,r3/2), 2: mp.polylog(3,r3/3)}
GD = [('3-s', 3-r3), ('s-1', r3-1), ('2s-3', 2*r3-3), ('s-2', r3-2), ('(s-1)/2', (r3-1)/2),
      ('(2s-3)/3', (2*r3-3)/3), ('2-s', 2-r3), ('(3-s)/2', (3-r3)/2), ('3-2s', 3-2*r3),
      ('4-2s', 4-2*r3), ('2-2s/3', 2-2*r3/3), ('2s-2', 2*r3-2)]
GDd = dict(GD)
Li2R = lambda z: mp.re(mp.polylog(2, mp.mpc(z)))
Li3R = lambda z: mp.re(mp.polylog(3, mp.mpc(z)))
def Creal(a,b):
    return mp.quad(lambda t: mp.log(abs(1-a*t))*mp.log(abs(1-b*t))/t, [0,1])
def SwpR(a,b):
    return Creal(a,b) + mp.log(abs(1-a))*(Li2R(b) - Li2R(b/a)) - Li3R(b)
gam = {}; dl = {}
for i in IDX:
    for j in IDX:
        if i==j: continue
        gam[(i,j)] = X[j]/(X[j]-X[i]); dl[(i,j)] = gam[(i,j)]*(1-X[i])
T3a = {}
for i in IDX:
    for j in IDX:
        for k in IDX:
            if len({i,j,k})==3 and j<k:
                T3a[(i,j,k)] = mp.quad(lambda tt: mp.log(1-X[j]*tt)*mp.log(1-X[k]*tt)/(1-X[i]*tt), [0,1])
SwG = {}; SwD = {}
for (i,j,k) in sorted(T3a):
    SwG[(i,j,k)] = SwpR(gam[(i,j)], gam[(i,k)])
    SwD[(i,j,k)] = SwpR(dl[(i,j)], dl[(i,k)])
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
for k,vv in GD: add('L3g_'+k, Li3R(vv))
reps = ['3-2s','s-2','(2s-3)/3','2-s','(s-1)/2','2s-3']
for k in reps:
    for nm,lv in [('l2',l2),('l3',l3),('lu',lu)]:
        add('L2g%s*%s'%(k,nm), Li2R(GDd[k])*lv)
print("B54 built, t=%ds"%(time.time()-t0))
# quotient B54 via LLL at 150 digits
N = len(vals); C = 10**120
rows = []
for i in range(N):
    rr_ = [0]*N + [int(mp.nint(C*vals[i]))]; rr_[i] = 1
    rows.append(rr_)
BM = sp.Matrix(rows).lll(delta=sp.Rational(3,4))
rels = []
for i in range(BM.rows):
    vc = [int(BM[i,j]) for j in range(N)]
    nrm = max(abs(x) for x in vc)
    if nrm==0 or nrm > 10**6: continue
    resid = abs(sum(mp.mpf(vc[q])*vals[q] for q in range(N)))
    if resid < mp.mpf('1e-110'): rels.append(vc)
print("B54 internal rels at 150d:", len(rels), "t=%ds"%(time.time()-t0))
Rm = sp.Matrix(rels)
Rr, piv = Rm.rref()
pivset=set(piv); keep=[i for i in range(N) if i not in pivset]
Bn = [names[i] for i in keep]; Bv = [vals[i] for i in keep]
print("B54_ind size:", len(keep), "pivots:", sorted(pivset))
# PSLQ each R
for (i,j,k) in sorted(T3a):
    R = T3a[(i,j,k)] - (SwG[(i,j,k)]-SwD[(i,j,k)])/X[i]
    rel = mp.pslq([R]+Bv, tol=mp.mpf('1e-130'), maxcoeff=10**8, maxsteps=20000)
    if rel is None:
        print("R%d_%d%d: NO REL"%(i,j,k))
    else:
        resid = rel[0]*R + sum(rel[q+1]*Bv[q] for q in range(len(Bv)))
        print("R%d_%d%d: resid=%s tcoef=%s maxc=%s sup=%s"%(i,j,k, mp.nstr(abs(resid),3),
              rel[0], max(abs(x) for x in rel), [(Bn[q],rel[q+1]) for q in range(len(Bv)) if rel[q+1]!=0][:14]))
print("elapsed", time.time()-t0)
