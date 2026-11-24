# step_70: verify OCT9c edge-0 presentation symbol vs OCT9d TW3[0]
import os, sys, pickle, itertools
from fractions import Fraction as Fr
base = "<project>"
sys.path.insert(0, os.path.join(base,"research/08_theory/OCT8_symbol_test"))
from oct8 import Quotient, S_Li3, vec_of_rational as vf
import sympy as sp
O9B = os.path.join(base,"research/08_theory/OCT9b_octagon_debug")
O9D = os.path.join(base,"research/08_theory/OCT9d_octagon_short")
s0 = pickle.load(open(os.path.join(O9B,"sol_B0.pkl"),"rb"))
sf = pickle.load(open(os.path.join(O9D,"short_form.pkl"),"rb"))
TW3 = pickle.load(open(os.path.join(O9D,"TW3.pkl"),"rb"))

def S_mLi2logp(w, p):
    out = {}
    for a, ca in vf(Fr(1)-Fr(w)).items():
        for b, cb in vf(Fr(w)).items():
            for c, cc in vf(Fr(p)).items():
                coef = ca*cb*cc
                for t in [(a,b,c),(a,c,b),(c,a,b)]:
                    out[t] = out.get(t, Fr(0)) + coef
    return {k:v for k,v in out.items() if v}
def S_logprod(primes):
    opts = [[(l, c) for l, c in vf(Fr(q)).items()] for q in primes]
    out = {}
    for combo in itertools.product(*opts):
        letters = tuple(l for l, c in combo)
        coef = Fr(1)
        for l, c in combo: coef *= c
        for perm in set(itertools.permutations(letters)):
            out[perm] = out.get(perm, Fr(0)) + coef
    return {k:v for k,v in out.items() if v}

Q = Quotient((2,3,5,7))
def coord(S):
    v = sp.zeros(Q.dim3,1)
    for t,c in S.items():
        v[Q.bidx[t],0] += sp.Rational(c.numerator, c.denominator)
    return v

sols = [Fr(s) for s in s0['meta']['sols']]
solB0 = [Fr(s) for s in s0['solB0']]
negw = [Fr(s) for s in s0['meta']['negw']]
logprods = [tuple(lp) for lp in s0['meta']['logprods']]
nL = len(sols)
nP = len(negw)*4
print("nL =", nL, " products:", nP, " logprods:", len(logprods), " solB0 len:", len(solB0))
tot = sp.zeros(Q.dim3,1)
for c,z in zip(solB0[:nL], sols):
    if c: tot += sp.Rational(c.numerator,c.denominator)*coord(S_Li3(z, vf))
idx = nL
for w in negw:
    for p in (2,3,5,7):
        c = solB0[idx]; idx += 1
        if c: tot += sp.Rational(c.numerator,c.denominator)*coord(S_mLi2logp(w,p))
for lp in logprods:
    c = solB0[idx]; idx += 1
    if c: tot += sp.Rational(c.numerator,c.denominator)*coord(S_logprod(lp))
t0 = sp.Matrix(TW3[0]); t0 = t0.T if t0.shape[0]==1 else t0
diff = tot - t0
print("symbol of OCT9c edge-0 presentation minus TW3[0]: nonzero entries:", sum(1 for x in diff if x != 0))
MD = sp.Matrix.hstack(*[d.T for d in Q.D])
print("diff in D (decomposables):", sp.Matrix.hstack(MD, diff).rank() == MD.rank())
