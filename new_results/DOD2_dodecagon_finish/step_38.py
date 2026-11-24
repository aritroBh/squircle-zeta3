import mpmath as mp, pickle, time
import sympy as sp
mp.mp.dps = 150
t0=time.time()
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
N = len(vals); C = 10**120
rows = []
for i in range(N):
    rr_ = [0]*N + [int(mp.nint(C*vals[i]))]; rr_[i] = 1
    rows.append(rr_)
BM = sp.Matrix(rows).lll(delta=sp.Rational(3,4))
print("LLL done t=%ds"%(time.time()-t0))
nfound=0
for i in range(BM.rows):
    vc = [int(BM[i,j]) for j in range(N)]
    nrm = max(abs(x) for x in vc)
    resid = abs(sum(mp.mpf(vc[q])*vals[q] for q in range(N)))
    lc = abs(int(BM[i,N]))
    if i < 12:
        print("row %2d: maxc=%.1e resid=%.1e lastcol=%.1e" % (i, nrm, resid, lc))
    if nrm>0 and nrm <= 10**6 and resid < mp.mpf('1e-110'):
        nfound+=1
print("passing filter:", nfound)
