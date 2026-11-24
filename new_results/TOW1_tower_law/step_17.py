lib = '''
"""towlib: tower data at v=+i, depth-indexed, for the TOW1 job."""
import sys, os, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
WTS = os.path.normpath(os.path.join(HERE, '..', 'WTS_weight_k_structure', 'theory'))
from fractions import Fraction as Fr
from fastred import *

def q4(q): return (q.p0, q.p1, q.p2, q.p3)
def unq4(t): return QSI(*t)

_towers = {}
def tower(k):
    if k not in _towers:
        _towers[k] = pickle.load(open(os.path.join(WTS, f'tower_{k}.pkl'), 'rb'))
    return _towers[k]

def hcoef(i, j, k, r):
    """coefficient of (v-i)^{-r} in HigherPart(W_{ij}^{(k)}) at +i, as QSI; 0 if absent."""
    H = tower(k).get((i,j))
    if H is None: return ZERO
    return unq4(H.get(r, (0,0,0,0)))

def deg(key): return key[0]+key[1]

def u(i, j, k, delta):
    """depth-delta tower coefficient: order r = 2d-2-delta."""
    d = i+j
    r = 2*d-2-delta
    if r < 2: raise ValueError('delta too large')
    return hcoef(i, j, k, r)

_binteg = {}
def rho(key, delta):
    """initial data: depth-delta coefficient of the level-d integrand (HigherPart of R_ij)."""
    d = deg(key)
    if d not in _binteg:
        _binteg[d] = pickle.load(open(os.path.join(WTS, f'binteg_{d}.pkl'), 'rb'))
    P = _binteg[d][key]
    # P is a PF [real, ipos, ineg, const] of QSI (as loaded); higher part at +i: ipos[1:]
    ipos = P[1]
    r = 2*d-2-delta
    if r-1 < len(ipos):
        return ipos[r-1]
    return ZERO

def keys_of_degree(d):
    return [(i, d-i) for i in range(d, (d-1)//2, -1)]

# Taylor coefficients of DL1, DL2 at v=i: tau1[n] = DL1^{(n+1)}(i)/(n+1)!
# DL1 = 1/(v-b)+1/(v+a)-1/v  =>  DL1^{(s)}(i)/s! = (-1)^s[1/(i-b)^{s+1}+1/(i+a)^{s+1}-1/i^{s+1}]
def tau1(n):
    s = n+1
    return ((-1)**s) * ( (IQ-BQ)**(-s) + (IQ+AQ)**(-s) - IQ**(-s) )
def tau2(n):
    s = n+1
    return ((-1)**s) * ( (IQ-AQ)**(-s) + (IQ+BQ)**(-s) - IQ**(-s) )
'''
open('towlib.py','w').write(lib)
import importlib, sys
sys.path.insert(0,'.')
import towlib
from towlib import *
print('tau1_0, tau2_0 =', tau1(0), tau2(0))   # expect (-1-i)/2, (-1+i)/2
print('rho(2,0; d=0):', rho((2,0),0), ' expect -sqrt2 i/4')
print('rho(1,1; d=0):', rho((1,1),0), ' expect sqrt2 i/4')
print('u(2,0,k=4,d=0):', u(2,0,4,0))
print('keys deg 3:', keys_of_degree(3))
