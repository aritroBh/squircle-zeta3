import sys, os, pickle
from fractions import Fraction as Fr
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
os.nice(10)
import sympy as sp
data = pickle.load(open('e61_fit_system.pkl','rb'))
rows0, _ = data['rows'], data['mons']
# rows0 targets were scaled by Pi_3. Rebuild with generic denominator from raw: recompute targets.
# Simpler: reuse rows0 (target*Pi3). For denominator D(d): target' = target * D/Pi3 = target * D/Pi3.
def try_fit(deg, den):
    mons = [(i,j) for i in range(0,deg+1) for j in range(0,deg+1-i)]
    N = len(mons)
    A = []; b = []
    for key_row, (tre, tim) in rows0:
        # key_row = mons^old; I stored row as dd^i p^j with old mons (deg5). Reconstruct dd, p? Not stored.
        pass
    return None
# I need dd,p per row — rebuild the system from scratch with dd,p saved this time.
import pickle as pk
H = pk.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
from math import comb
def u(i,j,k,delta):
    return tuple(H[(i,j)][delta][k]) if (i,j) in H and delta in H[(i,j)] and k in H[(i,j)][delta] else (Fr(0),)*4
def qsi_add(x,y): return tuple(a+b for a,b in zip(x,y))
def qsi_scale(x,f): return tuple(a*f for a in x)
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gdiv(x,y):
    den2 = y[0]**2+y[1]**2
    return ((x[0]*y[0]+x[1]*y[1])/den2, (x[1]*y[0]-x[0]*y[1])/den2)
def qsi_div(x, y):
    A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
    C2m2D2 = (gmul(C,C)[0]-2*gmul(D,D)[0], gmul(C,C)[1]-2*gmul(D,D)[1])
    numA = (gmul(A,C)[0]-2*gmul(B,D)[0], gmul(A,C)[1]-2*gmul(B,D)[1])
    numB = (gmul(B,C)[0]-gmul(A,D)[0], gmul(B,C)[1]-gmul(A,D)[1])
    return (gdiv(numA, C2m2D2), gdiv(numB, C2m2D2))
def e_ratio(key, delta, s):
    i,j = key
    k0 = (2*i+1) if i==j else i+j
    out = (Fr(0),)*4
    for t in range(0, s+1):
        out = qsi_add(out, qsi_scale(u(i,j,k0+t,delta), Fr((-1)**(s-t)*comb(s,t))))
    return qsi_div(out, u(i,j,k0,0))
samples = []
for key in H:
    if key=='_meta': continue
    p,q = key; dd = p+q
    if p==q or 2*dd-4 < 6 or not (6 <= dd <= 30): continue
    samples.append((dd, p, e_ratio(key, 6, 1)[0]))
samples = samples[:80]
print('samples', len(samples))
pk.dump(samples, open('e61_samples.pkl','wb'))
