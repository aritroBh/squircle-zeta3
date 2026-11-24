import sys, os, pickle
from fractions import Fraction as Fr
from math import comb
base = '<project>/research/08_theory'
import sympy as sp
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
s2 = sp.sqrt(2); ii = sp.I
aq = 1+s2; bq = s2-1
_wc = {}
def w_parts(mm):
    if mm in _wc: return _wc[mm]
    w = sp.radsimp(sp.expand((-aq-ii)**(-mm) - (bq-ii)**(-mm)))
    wre = sp.expand(sp.re(w)); wim = sp.expand(sp.im(w))
    out = ((Fr(sp.Rational(wre.coeff(s2,0))), Fr(sp.Rational(wre.coeff(s2,1)))),
           (Fr(sp.Rational(wim.coeff(s2,0))), Fr(sp.Rational(wim.coeff(s2,1)))))
    _wc[mm] = out
    return out

# d=2 (2,0) top: M=2
key=(2,0); d=2; M=2
ikey, j = key
n = d+1
parent = (ikey+1, j)
e00 = u(parent[0], parent[1], n, 0)
E = ((e00[0], e00[2]), (e00[1], e00[3]))
print('E =', E)
Yr = (Fr(0),Fr(0)); Ys = (Fr(0),Fr(0))
for mm in range(1, 2*d):
    deltap = 2*n-2-(mm+1)
    if deltap < 0: break
    if M-1 > deltap//2: continue
    er = e_ratio(parent, deltap, M-1)
    print('mm',mm,'deltap',deltap,'er',er)
    erg = er[0]
    (wr0,wr1),(ws0,ws1) = w_parts(mm)
    print('  w parts:', (wr0,wr1),(ws0,ws1))
    f = Fr(-1, mm)
    Yr = (Yr[0] + f*(erg[0]*wr0 - erg[1]*wr1), Yr[1] + f*(erg[0]*wr1 + erg[1]*wr0))
    Ys = (Ys[0] + f*(erg[0]*ws0 - erg[1]*ws1), Ys[1] + f*(erg[0]*ws1 + erg[1]*ws0))
print('Yr,Ys =', Yr, Ys)
def ggmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
Zr = (ggmul(E[0],Yr)[0] + 2*ggmul(E[1],Ys)[0], ggmul(E[0],Yr)[1] + 2*ggmul(E[1],Ys)[1])
Zs = (ggmul(E[0],Ys)[0] + ggmul(E[1],Yr)[0], ggmul(E[0],Ys)[1] + ggmul(E[1],Yr)[1])
print('Zr,Zs =', Zr, Zs)
fac = -(ikey+1)
print('result:', fac*2*Zr[0], fac*2*Zs[0])
