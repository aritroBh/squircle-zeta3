import sys, os, time
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
from math import comb

p = sp.Symbol('p', integer=True, positive=True)
d = 2*p-1
i = sp.I
dd_sym, mm_sym = T.d, T.m

def eR(delta, s, dd, mm, conj=False):
    """e(delta,s)/e00 closed form (off-diagonal) at (dd,mm); conj: i -> -i."""
    ii = -sp.I if conj else sp.I
    if s == 0: val = T.arho(delta).subs({dd_sym: dd, mm_sym: mm, sp.I: ii})
    elif delta == 2 and s == 1: val = -1/(2*dd-3)
    elif delta == 3 and s == 1: val = ii*(dd-3)/(2*dd-3)
    elif delta == 4 and s == 1:
        val = (dd**3 - 9*dd**2 + sp.Rational(113,4)*dd - ii*dd/2 + ii*mm - 29)/((2*dd-3)*(2*dd-5))
    elif delta == 4 and s == 2: val = sp.Rational(2)/((2*dd-3)*(2*dd-5))
    elif delta == 5 and s == 2: val = -2*ii*(dd-4)/((2*dd-3)*(2*dd-5))
    else: raise ValueError((delta,s))
    return sp.expand(val)

def tau1(n, conj=False):
    if conj: return (-1)**n * T.tau2(n)
    return T.tau1(n)
def tau2(n, conj=False):
    if conj: return (-1)**n * T.tau1(n)
    return T.tau2(n)

rA0 = i*(4*p-3)/(p+1)
rB0 = (4*p-1)*(4*p-3)/((p+1)*p)
def rA0c(conj=False): return -rA0 if conj else rA0
def rB0c(conj=False): return rB0

def ehat(delta, s, conj=False):
    """ehat_pp(delta,s)/e00_{p,p-1}: parent (p+1,p) degree 2p+1."""
    tot = 0
    for dpp in range(0, delta+1):
        tot += tau1(delta-dpp, conj) * eR(dpp, s, 2*p+1, p+1, conj) * rB0c(conj) / (4*p-dpp-1)
    return sp.expand((p+1)*tot)

def Ediag(delta, s, conj=False):
    tot = 0
    t = 0
    while s+t <= delta//2:
        tot += (-1)**t * ehat(delta, s+t, conj)
        t += 1
    return sp.expand(tot)

def feeds(delta, s):
    """Lambda^(delta)(s) for adjacent key (p,p-1), normalized by e00."""
    tot = 0
    # DL1: off-diagonal parent (p+1,p-1), degree 2p
    for dp in range(0, delta+1):
        tot += (p+1)*tau1(delta-dp)*eR(dp, s, 2*p, p+1)*rA0/(2*d-dp-1)
    # DL2: diagonal parent (p,p) via Ediag conversion
    for dp in range(0, delta+1):
        tot += p*tau2(delta-dp)*Ediag(dp, s)/(2*d-dp-1)
    # fold
    mm_ = 4*p-4-delta
    for mu in range(mm_, 4*p-3):
        g = 0
        for dp in range(0, 4*p-4-mu+1):
            g -= Ediag(dp, s, conj=True)*tau1(4*p-4-mu-dp, conj=True)/(4*p-3-dp)
        tot += -p * sp.expand(g) * (-i)**mu * comb(mu-2, mu-mm_) * i**(mm_-2)
    return sp.expand(tot)

identities = [(3,0),(3,1),(4,0),(4,1),(4,2)]
for (delta, s) in identities:
    lhs = eR(delta, s, d, p) + (eR(delta, s+1, d, p) if s+1 <= delta//2 else 0)
    rhs = feeds(delta, s)
    t0 = time.time()
    rem = sp.cancel(sp.together(lhs - rhs))
    print(f'adjacent ({delta},{s}) remainder: {rem}   ({round(time.time()-t0,1)}s)', flush=True)
