common = r'''import mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre

_glq = {}
def gl_nodes(level, prec):
    """Gauss-Legendre rule with 3*2**(level-1) nodes on [-1,1]."""
    key = (level, prec)
    if key not in _glq:
        old = mp.mp.prec
        mp.mp.prec = prec
        q = GaussLegendre(mp.mp)
        pairs = q.calc_nodes(level, prec)
        xs = [p[0] for p in pairs]; ws = [p[1] for p in pairs]
        _glq[key] = (xs, ws)
        mp.mp.prec = old
    return _glq[key]

def glint(f, segs, level, prec=None):
    if prec is None: prec = mp.mp.prec
    x, w = gl_nodes(level, prec)
    tot = mp.mpf(0)
    for a, b in segs:
        a = mp.mpf(a); b = mp.mpf(b)
        c = (a+b)/2; h = (b-a)/2
        s = mp.mpf(0)
        for xi, wi in zip(x, w):
            s += wi*f(c + h*xi)
        tot += h*s
    return tot

w4 = mp.expjpi(mp.mpf(1)/4)   # omega = e^{i pi/4}
aa = mp.pi/(2*mp.sqrt(2))     # a = pi/(2 sqrt(2))
def H(s):
    return mp.sqrt(mp.pi)/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf(3)/2)
def K(s):
    return 1/(2*mp.cos(mp.pi*s/2))
def E(s):
    return mp.expj(mp.pi*s/2)
def Qk(k, s):
    return s**k/(1+s)**(k+1)
def betaL(L):
    return mp.log(1+mp.exp(-L))
def _usegs(k, res):
    umax = max(mp.mpf(40), (k+30)/(res+1)*3)
    return [[0,1],[1,3],[3,8],[8,20],[20,umax]]
def Pk(k, s, level=5):
    def f(u):
        return mp.exp(-(s+1)*u)*(betaL(u)-s*u)**k/(1+mp.exp(-u))
    return (-1)**k/mp.factorial(k)*glint(f, _usegs(k, mp.re(s)), level)
def Ak(k, s, level=5):
    def f(L):
        return mp.exp(-(s+1)*L)*((s*L-betaL(L))**k/(1+mp.exp(-L)) - (s*L)**k)
    return glint(f, _usegs(k, mp.re(s)), level)/mp.factorial(k)
def Jk(k):
    return (-1)**k*mp.log(2)**(k+1)/mp.factorial(k+1)
import re
def load_table():
    tab = {}
    for line in open('<project>/research/08_theory/D13f_sweep/mk_deg9.tsv'):
        p = line.rstrip('\n').split('\t')
        if p[0].isdigit():
            tab[int(p[0])] = mp.mpf(re.match(r'\[([-0-9.eE+]+)', p[1]).group(1))
    return tab
GAMMA_SEGS = [[0,0.25],[0.25,0.5],[0.5,0.75],[0.75,1.0],[1.0,1.5],[1.5,2.0],[2.0,3.0],
              [3.0,5.0],[5.0,8.0],[8.0,12.0],[12.0,20.0],[20.0,32.0],[32.0,52.0],[52.0,80.0]]
RAY_SEGS = list(GAMMA_SEGS)
'''
open('cx6r_common.py','w').write(common)
import importlib, mpmath as mp
import cx6r_common as C
importlib.reload(C)
mp.mp.dps = 60
v1 = C.glint(lambda x: x**35, [[0,2]], 4)
v2 = C.glint(lambda x: mp.exp(-3*x), [[0,1],[1,5]], 6)
print("GL x^35 err:", mp.nstr(abs(v1 - mp.mpf(2)**36/36), 5))
print("GL exp err:", mp.nstr(abs(v2 - (1-mp.exp(-15))/3), 5))
import time
t=time.time()
p = C.Pk(3, mp.mpf('0.5')+2*C.w4)
print("Pk(3, 0.5+2w) =", mp.nstr(p, 30), f"({time.time()-t:.2f}s)")
