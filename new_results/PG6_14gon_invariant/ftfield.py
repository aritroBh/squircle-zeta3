
"""Exact arithmetic in F = Q(a), a = 2 cos(2 pi/7), minpoly a^3 + a^2 - 2 a - 1 = 0.
Element: tuple (A,B,C,d) with A,B,C Fractions, d positive int: value (A + B a + C a^2)/d.
Reduced: gcd(A,B,C,d)=1 (on numerators), d>0.
Embeddings: a1 = 2cos(2pi/7) ~ 1.24698, a2 = 2cos(4pi/7) ~ -0.44504, a3 = 2cos(6pi/7) ~ -1.80194.
"""
from fractions import Fraction as Fr
from math import gcd

def FT(A=0, B=0, C=0, d=1):
    A,B,C,d = Fr(A),Fr(B),Fr(C),int(d)
    if d<0: A,B,C,d = -A,-B,-C,-d
    g = 1
    for x in (A.numerator, B.numerator, C.numerator):
        g = gcd(g, abs(x))
    den_lcm = A.denominator*B.denominator*C.denominator  # cheap lcm upper bound
    # normalize: make denominators common
    return (A,B,C,d)

def ft(A=0,B=0,C=0,d=1):
    return (Fr(A),Fr(B),Fr(C),int(d))

ZERO=(Fr(0),Fr(0),Fr(0),1); ONE=(Fr(1),Fr(0),Fr(0),1); AL=(Fr(0),Fr(1),Fr(0),1)

def fadd(x,y):
    A1,B1,C1,d1=x; A2,B2,C2,d2=y
    from math import lcm
    d = lcm(d1,d2)
    return (A1*(d//d1)+A2*(d//d2), B1*(d//d1)+B2*(d//d2), C1*(d//d1)+C2*(d//d2), d)
def fneg(x): A,B,C,d=x; return (-A,-B,-C,d)
def fsub(x,y): return fadd(x,fneg(y))
def fmul(x,y):
    A1,B1,C1,d1=x; A2,B2,C2,d2=y
    # a^3 = -a^2+2a+1 ; a^4 = 3a^2+a-1
    A = A1*A2 + (B1*C2+C1*B2) - C1*C2
    B = A1*B2+B1*A2 + 2*(B1*C2+C1*B2) - C1*C2
    C = A1*C2+C1*A2 + B1*B2 - (B1*C2+C1*B2) + 3*C1*C2
    return (A,B,C,d1*d2)
def fscale(x,q):
    A,B,C,d=x; q=Fr(q)
    return (A*q.numerator, B*q.numerator, C*q.numerator, d*q.denominator)
def finv(x):
    A,B,C,d = x
    # multiplication matrix of (A+B a+C a^2), columns = images of 1,a,a^2:
    # *1 = A + B a + C a^2 ; *a = C + (A+2C) a + (B-C) a^2 ; *a^2 = (B-C) + (2B-C) a + (A-B+3C) a^2
    import sympy as sp
    M = sp.Matrix([[sp.Rational(A), sp.Rational(C), sp.Rational(B-C)],
                   [sp.Rational(B), sp.Rational(A+2*C), sp.Rational(2*B-C)],
                   [sp.Rational(C), sp.Rational(B-C), sp.Rational(A-B+3*C)]])
    Minv = M.inv()
    col = Minv * sp.Matrix([sp.Rational(d),0,0])
    return (Fr(col[0]), Fr(col[1]), Fr(col[2]), 1)
def fdiv(x,y): return fmul(x, finv(y))
def fpow(x,k):
    if k<0: return fpow(finv(x),-k)
    r=ONE
    for _ in range(k): r=fmul(r,x)
    return r
def feq(x,y): return fsub(x,y)[0]==0 and fsub(x,y)[1]==0 and fsub(x,y)[2]==0
def fnorm(x):
    # resultant of (A+B t+C t^2)/d with t^3+t^2-2t-1
    import sympy as sp
    t = sp.symbols('t')
    A,B,C,dd = x
    f = sp.Rational(A)+sp.Rational(B)*t+sp.Rational(C)*t**2
    g = t**3+t**2-2*t-1
    r = sp.resultant(f,g,t)
    return Fr(r) / (dd**3)
import mpmath as mp
_EMB = None
def get_EMB(prec=60):
    global _EMB
    old = mp.mp.dps; mp.mp.dps = prec
    _EMB = [2*mp.cos(2*mp.pi/7), 2*mp.cos(4*mp.pi/7), 2*mp.cos(6*mp.pi/7)]
    mp.mp.dps = old
    return _EMB
EMB = get_EMB()
def feval(x, prec=50):
    old = mp.mp.dps; mp.mp.dps = max(prec, 30)
    A,B,C,d = x
    a = EMB[0]
    v = (mp.mpf(int(A.numerator))/A.denominator + mp.mpf(int(B.numerator))/B.denominator*a + mp.mpf(int(C.numerator))/C.denominator*a*a)/d
    mp.mp.dps = old
    return v
def feval_emb(x, k):
    A,B,C,d = x
    a = EMB[k]
    return (mp.mpf(int(A.numerator))/A.denominator + mp.mpf(int(B.numerator))/B.denominator*a + mp.mpf(int(C.numerator))/C.denominator*a*a)/d
def from_pslq(y, maxcoeff=10**6, dps=None):
    """reconstruct (A,B,C,d) with y = (A+B a+C a^2)/d from mpf y."""
    old = mp.mp.dps
    if dps: mp.mp.dps = dps
    a = 2*mp.cos(2*mp.pi/7)
    rel = mp.pslq([y, 1, a, a*a], maxcoeff=maxcoeff)
    mp.mp.dps = old
    if rel is None: return None
    q0,q1,q2,q3 = rel  # q0 y + q1 + q2 a + q3 a^2 = 0
    if q0 == 0: return None
    return (Fr(-q1,q0), Fr(-q2,q0), Fr(-q3,q0), 1)

def norm_t(x):
    """reduce (A,B,C,d) to canonical form."""
    from math import lcm, gcd
    A,B,C,d = x
    L = lcm(lcm(A.denominator,B.denominator), lcm(C.denominator,d))
    a,b,c,l = A*L, B*L, C*L, L
    assert a.denominator==1 and b.denominator==1 and c.denominator==1
    a,b,c = a.numerator, b.numerator, c.numerator
    g = gcd(gcd(gcd(abs(a),abs(b)),abs(c)), l)
    a,b,c,l = a//g, b//g, c//g, l//g
    return (Fr(a,l), Fr(b,l), Fr(c,l), 1)
