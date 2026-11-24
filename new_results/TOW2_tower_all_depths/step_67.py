import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
from towlib import *
from math import comb as cb
from fractions import Fraction as Fr

p = sp.Symbol('p')
i = sp.I
dd_sym, mm_sym = T.d, T.m

def eR(delta, s, dd, mm, conj=False):
    if s > delta//2: return sp.Integer(0)
    ii = -sp.I if conj else sp.I
    if s == 0: val = T.arho_simple(delta).subs({dd_sym: dd, mm_sym: mm, sp.I: ii})
    elif delta == 2 and s == 1: val = -1/(2*dd-3)
    elif delta == 3 and s == 1: val = ii*(dd-3)/(2*dd-3)
    else: raise ValueError((delta,s))
    return val

def tau1(n, conj=False):
    return (-1)**n * T.tau2f(n) if conj else T.tau1f(n)
def tau2(n, conj=False):
    return (-1)**n * T.tau1f(n) if conj else T.tau2f(n)

rA0 = i*(4*p-3)/(p+1)
rB0 = (4*p-1)*(4*p-3)/((p+1)*p)

def ehat(delta, s, conj=False):
    tot = 0
    for dpp in range(0, delta+1):
        tot += tau1(delta-dpp, conj) * eR(dpp, s, 2*p+1, p+1, conj) * rB0 / (4*p-dpp-1)
    return sp.expand((p+1)*tot)

def Ediag(delta, s, conj=False):
    tot = 0; t = 0
    while s+t <= delta//2:
        tot += (-1)**t * ehat(delta, s+t, conj)
        t += 1
    return sp.expand(tot)

P = 3; dd_ = 2*P-1
def e_from_data(key, delta, s):
    i,j = key; dd = i+j
    base = (2*i+1) if i==j else dd
    return sum( u(i,j,base+t, delta).scale(Fr((-1)**(s-t)*cb(s,t))) for t in range(0,s+1) )
def Ediag_data(key, delta, sigma):
    i,j = key
    if i != j: return e_from_data(key, delta, sigma)
    out = ZERO; t = 0
    while sigma+t <= delta//2:
        out = out + e_from_data(key, delta, sigma+t).scale(Fr((-1)**t))
        t += 1
    return out
from towlib import tau1 as t1q, tau2 as t2q
e00 = e_from_data((P,P-1),0,0)
e00s = e00.sympy()

# sympy pieces (normalized), evaluated at p=3, times e00
def ff(x, c):
    r = sp.Integer(1)
    for u in range(c): r *= (x-u)
    return r
dl1_sp = sum( (P+1)*tau1(3-dp)*eR(dp, 1, 2*P, P+1)*(i*(4*P-3)/(P+1))/(2*dd_-dp-1) for dp in range(0,4) )
dl2_sp = sum( P*tau2(3-dp)*Ediag(dp, 1).subs(p,P)/(2*dd_-dp-1) for dp in range(0,4) )
fold_sp = 0
for t in range(0, 4):
    mu = 4*p-4-t
    g = 0
    for dp in range(0, t+1):
        g -= Ediag(dp, 1, conj=True)*tau1(t-dp, conj=True)/(4*p-3-dp)
    binom = ff(mu-2, 3-t)/sp.factorial(3-t)
    fold_sp += (-p * sp.expand(g) * (-i)**mu * binom * i**(4*p-4-3-2)).subs(p,P)
lhs_sp = eR(3,1,dd_,P)

# data pieces (absolute)
dl1_data = sum( e_from_data((P+1,P-1), dp, 1) * t1q(3-dp).scale(Fr(P+1, 2*dd_-dp-1)) for dp in range(0,4) )
dl2_data = sum( Ediag_data((P,P), dp, 1) * t2q(3-dp).scale(Fr(P, 2*dd_-dp-1)) for dp in range(0,4) )
fold_data = ZERO
mm_ = 4*P-4-3
for mu in range(mm_, 4*P-3):
    g = ZERO
    for dp in range(0, 4*P-4-mu+1):
        g = g - Ediag_data((P,P), dp, 1) * t1q(4*P-4-mu-dp).scale(Fr(1, 4*P-3-dp))
    fold_data = fold_data + g.conj_i() * (IQ**(-mu)) * Fr(cb(mu-2, mu-mm_)) * (IQ**(mm_-2))
fold_data = fold_data.scale(-P)
lhs_data = e_from_data((P,P-1),3,1)

def show(name, spval, dataval):
    spabs = sp.simplify(spval*e00s) if isinstance(spval, sp.Expr) else spval
    print(f'{name}: sympy*e00 = {spabs}   data = {dataval.sympy()}')
show('DL1', dl1_sp, dl1_data)
show('DL2', dl2_sp, dl2_data)
show('FOLD', fold_sp, fold_data)
show('LHS e(3,1)', lhs_sp, lhs_data)
