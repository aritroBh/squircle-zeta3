import sys, os, pickle, json
from fractions import Fraction as Fr
from math import comb
base = '<project>/research/08_theory'
work = f'{base}/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow3lib as T3
H = pickle.load(open(f'{base}/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
WTS = {k: pickle.load(open(f'{base}/WTS_weight_k_structure/theory/tower_{k}.pkl','rb')) for k in range(2,15)}
def u(i,j,k,delta):
    if (i,j) in H and delta in H[(i,j)] and k in H[(i,j)][delta]:
        return tuple(H[(i,j)][delta][k])
    n = i+j; r = 2*n-2-delta
    rec = WTS.get(k, {}).get((i,j), {}).get(r)
    return tuple(rec) if rec else (Fr(0),)*4
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

dsym, msym = T3.d, T3.m
def e_closed(delta, s, dd, mm):
    v = T3.e(delta, s, dd, mm)
    return (Fr(sp.Rational(sp.re(v))), Fr(sp.Rational(sp.im(v))))
def e83(dd, mm):
    num = 20*dd**3 - 300*dd**2 + (1549-10*sp.I)*dd + 20*sp.I*mm - 2592
    den = 4*T3.Pi(4, dd)
    v = sp.expand(num/den)
    return (Fr(sp.Rational(sp.re(v))), Fr(sp.Rational(sp.im(v))))
def S_d(d):
    if d == 2: return e_closed(2, 0, 3, 3)
    if d == 3: return e_closed(4, 1, 4, 4)
    if d == 4: return e_closed(6, 2, 5, 5)
    if d == 5: return e83(6, 6)
    if d == 6: return e_ratio((7,0), 10, 4)[0]
    if d == 7: return e_ratio((8,0), 12, 5)[0]

A = T3.A
def Pi_val(l, dd):
    out = 1
    for jj in range(1, l+1): out *= (2*dd-2*jj-1)
    return out
def Md_compact(d):
    n = d+1
    S = S_d(d); s0, s1 = S
    o1 = Fr((-1)**(d-1)*A[d-1], Pi_val(d-2, n))
    e0 = Fr((-1)**(d-2)*A[d-1], Pi_val(d-2, n))
    if d % 2 == 0:
        bracket = Fr(s0-s1, 2) - Fr(o1, 4) - Fr(e0, 6)
    else:
        bracket = Fr(s0+s1, 2) - Fr(o1, 4) - Fr(e0, 12)
    return 4*((-1)**(d-1))*Pi_val(d-1, n)*bracket
expected = {2:-1, 3:10, 4:15, 5:131, 6:394, 7:3151}
print('d : S_d (s0, s1) : M_d compact vs expected')
for d in range(2, 8):
    S = S_d(d)
    Md = Md_compact(d)
    print(f'd={d}: S=({S[0]}, {S[1]})   M={Md}  exp={expected[d]}  {"OK" if Md==expected[d] else "FAIL"}')
