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
    # TOWN1 first, else WTS towers
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

# g_m = w_m / sqrt2 as Gaussian rational
s2 = sp.sqrt(2); ii = sp.I
aq = 1+s2; bq = s2-1
_gc = {}
def g_m(mm):
    if mm in _gc: return _gc[mm]
    w = sp.radsimp(sp.expand((-aq-ii)**(-mm) - (bq-ii)**(-mm)))
    wim = sp.expand(w/s2)
    re = sp.Rational(sp.re(wim)); im = sp.Rational(sp.im(wim))
    _gc[mm] = (Fr(re), Fr(im))
    return _gc[mm]

dsym, msym = T3.d, T3.m
def e_closed(delta, s, dd, mm):
    v = T3.e(delta, s, dd, mm)
    return (sp.Rational(sp.re(v)), sp.Rational(sp.im(v)))

def e83(dd, mm):
    num = 20*dd**3 - 300*dd**2 + (1549-10*sp.I)*dd + 20*sp.I*mm - 2592
    den = 4*T3.Pi(4, dd)
    v = sp.expand(num/den)
    return (sp.Rational(sp.re(v)), sp.Rational(sp.im(v)))

def subtop_value(d):
    """e(2d-2, d-2) at (d+1, d+1) as Gaussian rational."""
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

def Md_map(d):
    n = d+1
    e00 = u(d+1, 0, n, 0)
    Es = (e00[1], e00[3])  # sqrt2 part, Gaussian
    g1, g2, g3 = g_m(1), g_m(2), g_m(3)
    # tops at (n, n)
    etop_even = (Fr((-1)**(d-1) * A[d]), Fr(0)); etop_even = (etop_even[0]/Pi_val(d-1, n), Fr(0))
    odd_top = ((-1)**(d-1) * A[d-1])  # coefficient with i*(1): e_top(2d-3) at (n,n): (-1)^{d-1} i A_{d-1} (n-(d-2)-2)/Pi = (-1)^{d-1} i A_{d-1} (1)/Pi_{d-2}(n)
    etop_odd = (Fr(0), Fr((-1)**(d-1) * A[d-1], Pi_val(d-2, n)))
    etop_even2 = (Fr((-1)**(d-2) * A[d-1], Pi_val(d-2, n)), Fr(0))
    S = subtop_value(d)
    S = (Fr(S[0]), Fr(S[1]))
    # G_top = -etop_even * g1 ; G_sub = -S*g1 - (1/2) etop_odd*g2 - (1/3) etop_even2*g3
    def ggmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
    Gt = ggmul((-etop_even[0], -etop_even[1]), g1)
    Gs = ggmul((-S[0], -S[1]), g1)
    t2 = ggmul(etop_odd, g2); Gs = (Gs[0] - t2[0]/2, Gs[1] - t2[1]/2)
    t3 = ggmul(etop_even2, g3); Gs = (Gs[0] - t3[0]/3, Gs[1] - t3[1]/3)
    # R = Re[Es Gs]/Re[Es Gt]
    num = ggmul(Es, Gs)[0]; den = ggmul(Es, Gt)[0]
    R = num/den
    return 2*A[d]*R, R

expected = {2:-1, 3:10, 4:15, 5:131, 6:394, 7:3151, 8:13849, 9:119969, 10:688074}
for d in range(2, 8):
    Md, R = Md_map(d)
    print(f'd={d}: R_d(0) = {R} = {sp.nsimplify(R)},  M_d = {Md}   expected {expected[d]}   {"OK" if Md==expected[d] else "FAIL"}')
