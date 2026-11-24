"""LEAD2 shared library: E(y) coefficients, zigzag, kernel D-table, master integrals."""
import os, pickle
import mpmath
from mpmath import mpf, sqrt, log, pi, zeta, gamma, sin, cosh, sinh, exp
from fractions import Fraction as Fr
from math import comb

BASE = '<project>'
THEORY = os.path.join(BASE, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')

def E_coeffs(N):
    """Taylor coeffs e_0..e_N of E(y) = Gamma(1+y/2)Gamma(1/2-y/2)/sqrt(pi),
    via log E = y log2 + sum_{n>=2} s_n y^n/n, s_n = zeta(n) (n even), eta(n) (n odd)."""
    lam = log(2)
    L = [mpf(0)]*(N+2); L[1] = lam
    for n in range(2, N+2):
        L[n] = (zeta(n) if n % 2 == 0 else (1-mpf(2)**(1-n))*zeta(n))/n
    E = [mpf(0)]*(N+2); E[0] = mpf(1)
    for n in range(1, N+2):
        E[n] = sum(k*L[k]*E[n-k] for k in range(1, n+1))/n
    return E

def Efun(y):
    return gamma(1+y/2)*gamma(mpf('0.5')-y/2)/sqrt(pi)

def zigzag(N):
    """A_0..A_N: sec x + tan x = sum A_n x^n/n!."""
    import sympy as sp
    x = sp.symbols('x')
    ser = sp.series(sp.sec(x)+sp.tan(x), x, 0, N+2).removeO().expand()
    return [int(ser.coeff(x, n)*sp.factorial(n)) for n in range(N+2)]

def eps(m):
    return [1, -1, -1, 1][m % 4]

def load_kernels(Kmax=18):
    kern = {}
    for K in range(1, Kmax+1):
        with open(os.path.join(THEORY, f'kern_{K}.pkl'), 'rb') as f:
            kern[K] = pickle.load(f)
    return kern

POLE = {(0,0):"0",(1,1):"a",(-1,-1):"-a",(-1,1):"b",(1,-1):"-b"}

def extract_Dtop(kern, dmax=8):
    """Newton top differences over k = d+1..18 for every (key, pole). Exact Fractions."""
    def ckk(k, key, r):
        P = kern.get(k, {}).get(key)
        if P is None or r not in P['real']: return Fr(0)
        return P['real'][r][0]
    Dtop = {}
    for d in range(1, dmax+1):
        keys = sorted({key for K in kern for key in kern[K] if sum(key) == d}, key=lambda kk: -kk[0])
        for key in keys:
            poles = sorted({r for K in kern for r in kern[K].get(key, {'real': {}})['real']})
            rec = {}
            for r in poles:
                seq = [ckk(K, key, r) for K in range(d+1, 19)]
                D = []; s = list(seq)
                for m in range(1, d+1):
                    D.append(s[0]); s = [x2-x1 for x1, x2 in zip(s, s[1:])]
                rec[POLE[r]] = {"Dtop": D[d-1], "poly_ok": all(t == Fr(0) for t in s)}
            Dtop[(d, key)] = rec
    return Dtop

def delta_frac(d):
    """delta_d = [x^d](sech x + tanh x) as Fraction."""
    import sympy as sp
    x = sp.symbols('x')
    ser = sp.series((1+sp.sinh(x))/sp.cosh(x), x, 0, d+1).removeO()
    c = sp.Rational(ser.coeff(x, d))
    return Fr(int(c.p), int(c.q))

def predicted_D(d, i, j, pole, de):
    """The binomial x zigzag table (HYPOTHESIS for all d; exact for d<=9)."""
    if pole == '-a': return (-1)**(j+1)*comb(d, j)*de
    if d % 2 == 1:
        if pole == 'a': return (-1)**j*comb(d, j)*de*(Fr(1) if j < (d-1)//2 else Fr(1, 2))
        if j == (d-1)//2:
            if pole == '-b': return Fr((-1)**(j+1)*comb(d, j))*de/2
            if pole == '0': return Fr((-1)**j*comb(d, j))*de/2
    else:
        if pole == '-b' and j < d//2: return (-1)**j*comb(d, j)*de
        if pole == '0' and j == d//2: return Fr((-1)**j*comb(d, j))*de/2
    return Fr(0)

# ---- u-line measures: m_e, m_odd, m_oddlast; S = sqrt(1+u^2) ----
def m_e(u):
    S = sqrt(1+u**2)
    return (1 - sqrt(2)/S)/(1-u**2)
def m_odd(u):
    S = sqrt(1+u**2)
    return (1 + sqrt(2)*u/S)/(1-u**2)
def m_oddlast(u):
    return m_odd(u) - 1/(2*(1-u))
def L1(u): return log((1+u)/2)
def L2(u): return log((1-u)/2)
