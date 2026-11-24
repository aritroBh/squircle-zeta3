import sys, os, time, pickle
from fractions import Fraction as Fr
base = '<project>/research/08_theory'
work = base + '/TOW4_tower_depth8'
sys.path.insert(0, work)
import sympy as sp
dsym, msym = sp.symbols('d m')
H = pickle.load(open(base + '/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))

def u(i,j,k,delta):
    return tuple(H[(i,j)][delta][k]) if (i,j) in H and delta in H[(i,j)] and k in H[(i,j)][delta] else None
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gdiv(x,y):
    den = y[0]**2+y[1]**2
    return ((x[0]*y[0]+x[1]*y[1])/den, (x[1]*y[0]-x[0]*y[1])/den)
def qsi_div(x, y):
    A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
    C2m2D2 = (gmul(C,C)[0]-2*gmul(D,D)[0], gmul(C,C)[1]-2*gmul(D,D)[1])
    numA = (gmul(A,C)[0]-2*gmul(B,D)[0], gmul(A,C)[1]-2*gmul(B,D)[1])
    numB = (gmul(B,C)[0]-gmul(A,D)[0], gmul(B,C)[1]-gmul(A,D)[1])
    return (gdiv(numA, C2m2D2), gdiv(numB, C2m2D2))

acl = {7: sp.sympify(pickle.load(open(work+'/a7_closed.pkl','rb')), locals={'d':dsym,'m':msym,'I':sp.I}),
       8: sp.sympify(pickle.load(open(work+'/a8_closed.pkl','rb')), locals={'d':dsym,'m':msym,'I':sp.I})}

t0=time.time()
for delta in [7,8]:
    bad=0; n=0; skipped=0
    aa = acl[delta]
    for key in H:
        if key=='_meta': continue
        p,q = key; dd = p+q
        if 2*dd-4 < delta: continue
        k0 = (2*p+1) if p==q else dd
        num = u(p,q,k0,delta); den = u(p,q,k0,0)
        if num is None or den is None or den==(Fr(0),)*4:
            skipped+=1; continue
        rv = qsi_div(num, den)
        cv = sp.expand(aa.subs({dsym: dd, msym: p}))
        expected = (sp.Rational(sp.re(cv)), sp.Rational(sp.im(cv)))
        n += 1
        if rv[1] != (Fr(0),Fr(0)) or rv[0] != (Fr(expected[0]), Fr(expected[1])):
            bad += 1
            if bad<=3: print("BAD", key, rv, expected)
    print(f"rho({delta})/e00 == a_{delta}: tested {n} cells, bad {bad}, skipped(missing) {skipped}  ({round(time.time()-t0,1)}s)", flush=True)
