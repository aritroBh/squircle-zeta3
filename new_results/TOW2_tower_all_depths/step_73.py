import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)

script = '''"""prove_E3b34.py -- PROOF of diagonal blindness (E3) at depths delta = 3, 4, for all p.
The feed of Delta_pp(delta') = rho_pp(delta') - sum_s (-1)^s ehat_pp(delta',s)  (delta' <= delta)
into the child (p,p-1) at depth delta (DL2-direct + fold) is identically zero.
Exact rational identity in p over Q(i); remainder 0 = complete proof.
"""
import sympy as sp
import tow2lib as T

p = sp.Symbol('p')
i = sp.I
dd_sym, mm_sym = T.d, T.m

def eR(delta, s, dd, mm, conj=False):
    """e(delta,s)/e00 of the OFF-DIAGONAL parent (p+1,p) at (dd,mm)=(2p+1,p+1)."""
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

rB0d = 2*i*(4*p-1)/(p+1)          # e00_{p+1,p}/e00_{pp}  (diagonal ratio)
def rB0dc(conj=False): return -rB0d if conj else rB0d

def ehat(delta, s, conj=False):
    """ehat_pp(delta,s)/e00_{pp}."""
    tot = 0
    for dpp in range(0, delta+1):
        tot += tau1(delta-dpp, conj) * eR(dpp, s, 2*p+1, p+1, conj) * rB0dc(conj) / (4*p-dpp-1)
    return sp.expand((p+1)*tot)

def rho_pp(delta, conj=False):
    ii = -sp.I if conj else sp.I
    return T.arho_simple(delta).subs({dd_sym: 2*p, mm_sym: p, sp.I: ii})

def Delta(delta, conj=False):
    return sp.expand( rho_pp(delta, conj) - sum( (-1)**s * ehat(delta, s, conj) for s in range(0, delta//2+1) ) )

def Dc(delta):
    return -Delta(delta, conj=True)

def gfun(t):
    """conj(g_mu) for mu = 4p-4-t."""
    return sum( -Dc(dp)*tau1(t-dp, conj=True)/(4*p-3-dp) for dp in range(0, t+1) )

def feed(delta):
    direct = p*sum( tau2(delta-dp)*Delta(dp)/(4*p-3-dp) for dp in range(0, delta+1) )
    ft = 0
    for t in range(0, delta+1):
        # coefficient c_t = (-i)^{4p-4-t} * i^{4p-6-delta} * binom(4p-6-t, delta-t)
        c = (-i)**(-4-t) * i**(-6-delta)
        from math import comb
        binom = sp.binomial(4*p-6-t, delta-t) if False else None
        # symbolic binomial -> falling factorial
        f = sp.Integer(1)
        for u in range(delta-t): f *= (4*p-6-t-u)
        binom = f/sp.factorial(delta-t)
        ft += -p * sp.expand(gfun(t)) * c * binom
    return sp.expand(direct + ft)

import time
for delta in (3, 4):
    t0 = time.time()
    f = feed(delta)
    terms = [sp.cancel(t) for t in sp.Add.make_args(f)]
    dens = [sp.denom(t) for t in terms]
    DEN = sp.Integer(1)
    for dn in dens: DEN *= dn
    rem = sp.expand(f * DEN)
    is0 = sp.simplify(rem) == 0
    print(f'E3 blindness delta={delta}: remainder*DEN == 0 : {is0}  ({round(time.time()-t0,1)}s)', flush=True)
'''
open('prove_E3b34.py','w').write(script)
import subprocess, sys
r = subprocess.run(['nice','-n','10',sys.executable,'prove_E3b34.py'], capture_output=True, text=True, timeout=2800)
print(r.stdout, r.stderr[-300:] if r.returncode else '')
