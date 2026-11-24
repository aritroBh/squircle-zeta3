import mpmath
mpmath.mp.dps = 80
from mpmath import mpf, sqrt, log, quad, pi, gamma, zeta
from math import comb
S = lambda u: sqrt(1+u**2)
def m_odd(u): return (1 + sqrt(2)*u/S(u))/(1-u**2)
def m_even(u): return (1 - sqrt(2)/S(u))/(1-u**2)
def m_diag(u): return (1 - sqrt(2)/S(u))/(2*(1+u))
def m_oddlast(u): return 1/(2*(1+u)) + sqrt(2)*u/(S(u)*(1-u**2))
def Mcombo(d, j, i=None):
    if i is None: i = d-j
    if d%2==1:
        m = m_odd if j < (d-1)//2 else m_oddlast
    else:
        m = m_even if j < d//2 else m_diag
    return quad(lambda u: log((1+u)/2)**i * log((1-u)/2)**j * m(u), [-1,0,1])
# e_w
NN=12
lam = log(2)
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
L=[mpf(0)]*(NN+1); L[1]=lam
for n in range(2,NN+1): L[n]=(zeta(n) if n%2==0 else eta(n))/n
E=[mpf(0)]*(NN+1); E[0]=mpf(1)
for n in range(1,NN+1): E[n]=sum(k*L[k]*E[n-k] for k in range(1,n+1))/n
import sympy as sp
X = sp.symbols('X'); ser = sp.series((1+sp.sinh(X))/sp.cosh(X),X,0,10).removeO()
delta = {d: mpf(str(sp.Rational(ser.coeff(X,d)))) for d in range(1,9)}
for d in range(1,9):
    tot = sum((-1)**(j+1)*comb(d,j)*Mcombo(d,j) for j in range(0, d//2+1 if d%2==0 else (d-1)//2+1))
    w = d+1
    target = mpmath.factorial(d)/2 * E[w]
    print(f"d={d} (w={w}): total={mpmath.nstr(tot,25)}  (d!/2)e_{{d+1}}={mpmath.nstr(target,25)}  diff={mpmath.nstr(tot-target,8)}")
