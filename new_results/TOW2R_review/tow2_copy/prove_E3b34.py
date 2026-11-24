"""prove_E3b34.py -- PROOF of diagonal blindness (E3) at depths delta = 3, 4, all p."""
import sympy as sp, time
import tow2lib as T

p = sp.Symbol('p')
i = sp.I
dd_sym, mm_sym = T.d, T.m

def eR(delta, s, dd, mm, conj=False):
    if s > delta//2: return sp.Integer(0)
    ii = -sp.I if conj else sp.I
    if s == 0: val = T.arho_simple(delta).subs({dd_sym: dd, mm_sym: mm, sp.I: ii})
    elif delta == 2 and s == 1: val = -1/(2*dd-3)
    elif delta == 3 and s == 1: val = ii*(dd-3)/(2*dd-3)
    elif delta == 4 and s == 1:
        val = (dd**3 - 9*dd**2 + sp.Rational(113,4)*dd - ii*dd/2 + ii*mm - 29)/((2*dd-3)*(2*dd-5))
    elif delta == 4 and s == 2: val = sp.Rational(2)/((2*dd-3)*(2*dd-5))
    else: raise ValueError((delta,s))
    return val

def tau1(n, conj=False):
    return (-1)**n * T.tau2f(n) if conj else T.tau1f(n)
def tau2(n, conj=False):
    return (-1)**n * T.tau1f(n) if conj else T.tau2f(n)

rB0d = 2*i*(4*p-1)/(p+1)

def ehat(delta, s, conj=False):
    tot = 0
    for dpp in range(0, delta+1):
        tot += tau1(delta-dpp, conj) * eR(dpp, s, 2*p+1, p+1, conj) * (rB0d if not conj else -rB0d) / (4*p-dpp-1)
    return sp.expand((p+1)*tot)

def rho_pp(delta, conj=False):
    ii = -sp.I if conj else sp.I
    return T.arho_simple(delta).subs({dd_sym: 2*p, mm_sym: p, sp.I: ii})

def Delta(delta, conj=False):
    D = rho_pp(delta, conj) - sum( (-1)**s * ehat(delta, s, conj) for s in range(0, delta//2+1) )
    return sp.cancel(D)

def Dc(delta):
    return sp.expand(-Delta(delta, conj=True))

def gfun(t):
    return sum( -Dc(dp)*tau1(t-dp, conj=True)/(4*p-3-dp) for dp in range(0, t+1) )

def feed(delta):
    direct = p*sum( tau2(delta-dp)*Delta(dp)/(4*p-3-dp) for dp in range(0, delta+1) )
    ft = 0
    for t in range(0, delta+1):
        c = (-i)**(-4-t) * i**(-6-delta)
        f = sp.Integer(1)
        for u in range(delta-t): f *= (4*p-6-t-u)
        binom = f/sp.factorial(delta-t)
        ft += -p * sp.expand(gfun(t)) * c * binom
    return sp.expand(direct + ft)

import sys
delta = int(sys.argv[1])
t0 = time.time()
f = feed(delta)
rem = sp.together(f)
num = sp.fraction(rem)[0]
num = sp.expand(num)
is0 = num == 0
print(f'E3 blindness delta={delta}: numerator == 0 : {is0}  ({round(time.time()-t0,1)}s)', flush=True)
assert is0
