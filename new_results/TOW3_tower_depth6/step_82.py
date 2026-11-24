import sys, os, pickle
from fractions import Fraction as Fr
from math import comb
base = '<project>/research/08_theory'
work = f'{base}/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
H = pickle.load(open(f'{base}/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
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
    for (ii2,jj),(re,im) in zip(f['mons'], f['sol']):
        num += (sp.Rational(re) + sp.Rational(im)*sp.I) * dsym**ii2 * msym**jj
    return num / T.Pi(f['b'])
e83 = fit_to_sym(efits[(8,3)])
print('e(8,3) =', sp.cancel(e83))
bad = 0; n = 0
for key in H:
    if key=='_meta': continue
    p,q = key; dd = p+q
    if p==q or 2*dd-4 < 8 or dd > 40: continue
    rv = e_ratio(key, 8, 3)
    cv = sp.expand(e83.subs({dsym: dd, msym: p}))
    expected = (sp.Rational(sp.re(cv)), sp.Rational(sp.im(cv)))
    n += 1
    if rv[1] != (Fr(0),Fr(0)) or rv[0] != expected: bad += 1
print(f'e(8,3) fit vs data: {n} keys, bad {bad}')
