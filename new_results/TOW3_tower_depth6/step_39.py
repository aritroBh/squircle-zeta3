import sys, os, time, pickle
from fractions import Fraction as Fr
from math import comb
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
os.nice(10)
import sympy as sp
import tow2lib as T
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def u(i,j,k,delta):
    return tuple(H[(i,j)][delta][k]) if (i,j) in H and delta in H[(i,j)] and k in H[(i,j)][delta] else (Fr(0),)*4
def qsi_add(x,y): return tuple(a+b for a,b in zip(x,y))
def qsi_scale(x,f): return tuple(a*f for a in x)
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
def e_ratio(key, delta, s):
    i,j = key
    k0 = (2*i+1) if i==j else i+j
    out = (Fr(0),)*4
    for t in range(0, s+1):
        out = qsi_add(out, qsi_scale(u(i,j,k0+t,delta), Fr((-1)**(s-t)*comb(s,t))))
    return qsi_div(out, u(i,j,k0,0))

dsym, msym = T.d, T.m
efits = pickle.load(open('e_fits.pkl','rb'))
def fit_to_sym(f):
    num = sp.Integer(0)
    for (ii,jj),(re,im) in zip(f['mons'], f['sol']):
        num += (sp.Rational(re) + sp.Rational(im)*sp.I) * dsym**ii * msym**jj
    return num / T.Pi(f['b'])

# verify (5,1),(6,2) fits on ALL keys d<=40
for (delta,s) in [(5,1),(6,2)]:
    expr = fit_to_sym(efits[(delta,s)])
    bad = 0; n = 0
    for key in H:
        if key=='_meta': continue
        p,q = key; dd = p+q
        if p==q or 2*dd-4 < delta or dd > 40: continue
        rv = e_ratio(key, delta, s)
        cv = sp.expand(expr.subs({dsym: dd, msym: p}))
        expected = (sp.Rational(sp.re(cv)), sp.Rational(sp.im(cv)))
        n += 1
        if rv[1] != (Fr(0),Fr(0)) or rv[0] != expected: bad += 1
    print(f'fit ({delta},{s}) vs data: {n} keys, bad {bad}', flush=True)

# fit e(6,1): num = sum_{i+j<=5} c_ij d^i m^j, denominator Pi_3
mons = [(i,j) for i in range(0,6) for j in range(0,6-i)]
keys_fit = []
for key in H:
    if key=='_meta': continue
    p,q = key; dd = p+q
    if p==q or 2*dd-4 < 6: continue
    if 6 <= dd <= 30:
        keys_fit.append(key)
keys_fit = keys_fit[:80]
rows = []
for key in keys_fit:
    p,q = key; dd = p+q
    rv = e_ratio(key, 6, 1)
    assert rv[1] == (Fr(0),Fr(0)), ('sqrt2 part!', key)
    val = rv[0]  # Gaussian rational (re, im)
    pi3 = Fr((2*dd-3)*(2*dd-5)*(2*dd-7))
    target = (val[0]*pi3, val[1]*pi3)
    row = []
    for (i,j) in mons:
        row.append(dd**i * p**j)
    rows.append((row, target))
print('fit system:', len(rows), 'rows x', len(mons), 'mons')
pickle.dump({'rows': rows, 'mons': mons}, open('e61_fit_system.pkl','wb'))
