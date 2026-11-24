code = r'''
import mpmath as mp, pickle, time
import sympy as sp
mp.mp.dps = 620
t0 = time.time()
D1 = '../DOD1_dodecagon/'
r3 = mp.sqrt(3)
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
X = atoms['X']; T3a = atoms['T3a']
l2,l3,lu,z2,z3 = atoms['l2'],atoms['l3'],atoms['lu'],atoms['z2'],atoms['z3']
L2x,L3x,L3o = atoms['L2'],atoms['L3'],atoms['L3o']
g = pickle.load(open('gdswaps620.pkl','rb'))
SwG, SwD = g['SwG'], g['SwD']
GD = {'3-s': 3-r3, 's-1': r3-1, '2s-3': 2*r3-3, 's-2': r3-2, '(s-1)/2': (r3-1)/2,
      '(2s-3)/3': (2*r3-3)/3, '2-s': 2-r3, '(3-s)/2': (3-r3)/2, '3-2s': 3-2*r3,
      '4-2s': 4-2*r3, '2-2s/3': 2-2*r3/3, '2s-2': 2*r3-2}
L2g = {k: mp.re(mp.polylog(2, mp.mpc(vv))) for k,vv in GD.items()}
L3g = {k: mp.re(mp.polylog(3, mp.mpc(vv))) for k,vv in GD.items()}
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
for k in GD: add('L3g_'+k, L3g[k])
reps = ['3-2s','s-2','(2s-3)/3','2-s','(s-1)/2','2s-3']
for k in reps:
    for nm,lv in [('l2',l2),('l3',l3),('lu',lu)]: add('L2g%s*%s'%(k,nm), L2g[k]*lv)
NB = len(vals)
# the 12 T3a remainders
triples = sorted(T3a.keys())
for (i,j,k) in triples:
    R = T3a[(i,j,k)] - (SwG[(i,j,k)] - SwD[(i,j,k)])/X[i]
    add('R%d_%d%d'%(i,j,k), R)
N = len(vals)
print("basis:", NB, "total:", N, "t=%ds"%(time.time()-t0))
C = 10**540
rows = []
for i in range(N):
    rr_ = [0]*N + [int(mp.nint(C*vals[i]))]
    rr_[i] = 1
    rows.append(rr_)
M = sp.Matrix(rows)
print("LLL start", time.time()-t0)
B = M.lll(delta=sp.Rational(3,4))
print("LLL done", time.time()-t0)
found = {}
basis_rels = 0
for i in range(B.rows):
    vc = [int(B[i,j]) for j in range(N)]
    nrm = max(abs(x) for x in vc)
    if nrm==0 or nrm > 10**6: continue
    resid = abs(sum(mp.mpf(vc[q])*vals[q] for q in range(N)))
    if resid < mp.mpf('1e-500'):
        sup_R = [names[q] for q in range(NB,N) if vc[q]!=0]
        if sup_R:
            for nm in sup_R:
                found.setdefault(nm, []).append((vc, resid, nrm))
        else:
            basis_rels += 1
print("pure-basis relations:", basis_rels)
print("R-constants covered by a clean relation: %d/12"%len(found))
for nm in sorted(found):
    vc,resid,nrm = sorted(found[nm], key=lambda t:t[2])[0]
    sup = [(names[q],vc[q]) for q in range(N) if vc[q]!=0]
    print(nm, "maxc=%d resid=%s"%(nrm, mp.nstr(resid,3)))
    print("   ", sup)
pickle.dump({'names':names,'found':{k:[(vc,str(r),n) for vc,r,n in v] for k,v in found.items()}},
            open('t3a_remainder_rels.pkl','wb'))
print("elapsed", time.time()-t0)
'''
open('step_04.py','w').write(code)
exec(code)
