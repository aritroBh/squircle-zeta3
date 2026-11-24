# step_13 (fixed2): explicit CR-only Li3 presentations of each edge's TW3
import os, sys, itertools, pickle
from fractions import Fraction as Fr
base = "<project>"
sys.path.insert(0, os.path.join(base,"research/08_theory/OCT8_symbol_test"))
from oct8 import Quotient, S_Li3, vec_of_rational as vf
import sympy as sp

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

def toFr(x):
    r = sp.Rational(x); return Fr(int(r.p), int(r.q))

K12 = {0:(Fr(3,2),Fr(-7,12)), 1:(Fr(1,3),Fr(1,9)), 2:(Fr(5,3),Fr(-5,9)), 3:(Fr(1,4),Fr(1,24))}
TW3 = pickle.load(open(os.path.join(base,"research/08_theory/OCT9d_octagon_short/TW3.pkl"),"rb"))

CR_presentations = {}
for e, primes in [(0,(2,3,5,7)), (1,(2,3,5)), (2,(2,3,5)), (3,(2,3,5))]:
    k1,k2 = K12[e]; p1,p2 = -1/k1,-1/k2
    CR = cross_ratios_Q([Fr(0),Fr(1),None,p1,p2])
    Q = Quotient(primes)
    gens = [Q.coord(S_Li3(z, vf)) for z in CR]
    gens = [g.T if g.shape[0]==1 else g for g in gens]
    M = sp.Matrix.hstack(*gens + [d.T for d in Q.D])
    t = sp.Matrix(TW3[e]);  t = t.T if t.shape[0]==1 else t
    sol = M.gauss_jordan_solve(t)[0]
    ng = len(gens)
    coeffs = [toFr(sol[i,0]) for i in range(ng)]
    nz = [(c, z) for c,z in zip(coeffs, CR) if c != 0]
    Ssum = sp.zeros(Q.dim3,1)
    for c,z in nz:
        Ssum += sp.Rational(c.numerator,c.denominator)*gens[CR.index(z)]
    diff = t - Ssum
    MD = sp.Matrix.hstack(*[d.T for d in Q.D])
    inD = sp.Matrix.hstack(MD, diff).rank() == MD.rank()
    CR_presentations[e] = {'p1':p1,'p2':p2,'CR':CR,'presentation':nz,'diff_in_D':bool(inD)}
    print(f"edge {e}: {len(nz)} nonzero Li3 terms (of {len(CR)} CR points); TW3 - sum in D: {inD}")
    for c,z in nz: print(f"    {c} * Li3({z})")

pickle.dump(CR_presentations, open("CR_presentations.pkl","wb"))
print("saved CR_presentations.pkl")
