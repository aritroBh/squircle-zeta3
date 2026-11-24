# step_12: symbol-level span checks for the generic octagon edges
import os, sys, itertools, pickle
from fractions import Fraction as Fr
base = "<project>"
sys.path.insert(0, os.path.join(base,"research/08_theory/OCT8_symbol_test"))
import oct8
from oct8 import Quotient, S_Li3, vec_of_rational as vf
import sympy as sp

def S3orbit_q(z):
    z = Fr(z); pts=[z,1-z]
    if z!=0: pts.append(1/z)
    if z!=1: pts.append(1/(1-z))
    if z!=0 and z!=1: pts += [z/(z-1),(z-1)/z]
    return sorted(set(pts))

def cross_ratios_Q(points):
    out=set()
    for a,b,c,d in itertools.permutations(points,4):
        try:
            if a is None: val=(b-d)/(b-c)
            elif b is None: val=(a-c)/(a-d)
            elif c is None: val=(b-d)/(a-d)
            elif d is None: val=(a-c)/(b-c)
            else: val=(a-c)*(b-d)/((a-d)*(b-c))
        except ZeroDivisionError: continue
        out.add(Fr(val))
    return sorted(out)

K12 = {0:(Fr(3,2),Fr(-7,12)), 1:(Fr(1,3),Fr(1,9)), 2:(Fr(5,3),Fr(-5,9)), 3:(Fr(1,4),Fr(1,24))}
TW3 = pickle.load(open(os.path.join(base,"research/08_theory/OCT9d_octagon_short/TW3.pkl"),"rb"))
print("TW3 type:", type(TW3), "len", len(TW3))

for e, primes in [(1,(2,3,5)), (2,(2,3,5)), (3,(2,3,5)), (0,(2,3,5,7))]:
    k1,k2 = K12[e]; p1,p2 = -1/k1,-1/k2
    CR = cross_ratios_Q([Fr(0),Fr(1),None,p1,p2])
    Q = Quotient(primes)
    n = len(primes)
    print(f"=== edge {e}: letters {primes}; dim V^3 = {Q.dim3}, rank D = {Q.rank_D}, dim Q = {Q.dim3-Q.rank_D} (Witt (n^3-n)/3 = {(n**3-n)//3})")
    gens = [Q.coord(S_Li3(z, vf)) for z in CR]
    M = sp.Matrix.hstack(*[g.T for g in gens])
    # rank mod products: rank(gens + D) - rank(D)
    MD = sp.Matrix.hstack(M, *[d.T for d in Q.D])
    r_all = MD.rank(); rD = Q.rank_D
    print(f"  #CR = {len(CR)}, rank span S(Li3(CR)) mod D = {r_all-rD}")
    t = sp.Matrix(TW3[e]) if not isinstance(TW3[e], sp.Matrix) else TW3[e]
    if t.shape[0] == 1: t = t.T
    MDt = sp.Matrix.hstack(MD, t)
    print(f"  TW3[{e}] in span(S(Li3(CR_e)))+D : {MDt.rank()==r_all}")
