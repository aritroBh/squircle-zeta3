import mpmath as mp, pickle, time
import sympy as sp
mp.mp.dps = 620
t0=time.time()
D1 = '../DOD1_dodecagon/'
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
r3 = mp.sqrt(3)
X=atoms['X']; T3a=atoms['T3a']
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
NB=len(vals)
for (i,j,k) in sorted(T3a.keys()):
    add('R%d_%d%d'%(i,j,k), T3a[(i,j,k)] - (SwG[(i,j,k)]-SwD[(i,j,k)])/X[i])
N=len(vals)
# test pslq on ONE remainder to cross-check LLL
Rname = names[NB]
Rv = vals[NB]
print("pslq on", Rname, "t=%ds"%(time.time()-t0))
rel = mp.pslq([Rv]+vals[:NB], tol=mp.mpf('1e-580'), maxcoeff=10**8, maxsteps=20000)
print("pslq done t=%ds"%(time.time()-t0))
if rel is None:
    print("PSLQ: NO relation with maxcoeff 1e8")
else:
    resid = rel[0]*Rv + sum(rel[i+1]*vals[i] for i in range(NB))
    print("PSLQ relation, resid:", mp.nstr(abs(resid),4), "target coef:", rel[0])
    print([(names[i],rel[i+1]) for i in range(NB) if rel[i+1]!=0])
