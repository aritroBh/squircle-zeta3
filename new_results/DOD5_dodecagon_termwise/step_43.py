import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
mp.mp.dps = 620
cl = pickle.load(open('dod5_closures36.pkl','rb'))
ac = pickle.load(open(D2+"atom_coefs_pq.pkl","rb"))
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
from fractions import Fraction as Fr
# coefficient 1-parts and sqrt3-parts
c0map = {}; c1map = {}
for (kind, ij),(p,q) in ac.items():
    if kind=='Lam':
        c0map[('Lam',ij)] = mp.mpf(p.numerator)/p.denominator
        c1map[('Lam',ij)] = mp.mpf(q.numerator)/q.denominator
    else:
        i,j,k = ij
        xi_p, xi_q = psplit(Xs[i]); den = xi_p*xi_p - 3*xi_q*xi_q
        gp = (Fr(p)*xi_p - 3*Fr(q)*xi_q)/den; gq = (Fr(q)*xi_p - Fr(p)*xi_q)/den
        c0map[('SwG',ij)] = mp.mpf(gp.numerator)/gp.denominator; c1map[('SwG',ij)] = mp.mpf(gq.numerator)/gq.denominator
        c0map[('SwD',ij)] = -mp.mpf(gp.numerator)/gp.denominator; c1map[('SwD',ij)] = -mp.mpf(gq.numerator)/gq.denominator
rem0 = sum(c0map[k]*cl[k]['rem'] for k in cl)
rem1 = sum(c1map[k]*cl[k]['rem'] for k in cl)
print("sum c0*rem =", mp.nstr(rem0, 20), " (DOD4 rem0 was -579.1913644993620292824497)")
print("sum c1*rem =", mp.nstr(rem1, 20), " (DOD4 rem1 was 339.0352590365888847442647)")
# estimate the denominator size of the combined invisible relation
from math import lcm, gcd
dens = set()
for k in cl:
    rel = cl[k]['rel']
    for r in rel: dens.update([d for d in [abs(r)] if d>1])
print("distinct |coefs|>1 across all 36 relations:", len(dens))
import numpy as np
print("sample:", sorted(dens)[:20])
