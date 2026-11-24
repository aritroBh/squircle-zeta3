import sys, os, time
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
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
    return rho_pp(delta, conj) - sum( (-1)**s * ehat(delta, s, conj) for s in range(0, delta//2+1) )

t0 = time.time()
for delta in (3, 4):
    D = sp.cancel(Delta(delta))
    print(f'Delta({delta}): size {sp.count_ops(D)} ({round(time.time()-t0,1)}s)', flush=True)
    num, den = sp.fraction(D)
    print('   numerator factored:', sp.factor(num) if sp.count_ops(num) < 400 else f'big ({sp.count_ops(num)})', flush=True)
    print('   denominator:', sp.factor(den), flush=True)
