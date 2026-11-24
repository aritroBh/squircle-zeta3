
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, math, time, pickle
import mpmath as mp
import sympy as sp
import flint
import w5lib as W
import pipe
D=12
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
LL = pickle.load(open(pipe.cand_dir(D)+"/atom_lll.pkl","rb"))
names = LL["names"]; a = LL["exact_rels"][0][1]; a0 = a[0]
Z = LL["exact_zeros"]
N = len(cands); M = len(names) - 2 - N
c = [Fraction(-x, a0) for x in a]
Zf = [[Fraction(x) for x in z] for z in Z]
li5idx = list(range(2, 2+N))

def cost(cv): return max(abs(x.numerator) for x in cv[1:])
def nterms(cv): return sum(1 for x in cv[1:] if x)
def euclid2(cv): return sum(float(x*x) for x in cv[1:])

# LLL-reduce the zero lattice (integer rows)
Zi = [[int(x) for x in z] for z in Z]
R = flint.fmpz_mat(Zi).lll(gram='exact')
Zr = [[int(R[i,j]) for j in range(R.ncols())] for i in range(R.nrows())]
Zr = [z for z in Zr if any(z)]
print("reduced zero basis:", len(Zr), " max entry:", max(max(abs(x) for x in z) for z in Zr))
ZrF = [[Fraction(x) for x in z] for z in Zr]

cv = list(c)
print("start: terms", nterms(cv), "cost", cost(cv))

# phase 1: Babai nearest-plane (integer t), several passes with random shuffles
import random
def gso(B):
    n = len(B); m = len(B[0])
    Bs = [[Fraction(x) for x in row] for row in B]
    mu = [[Fraction(0)]*n for _ in range(n)]
    Bst = []
    norm = []
    for i in range(n):
        v = list(Bs[i])
        for j in range(i):
            d = sum(Bs[i][k]*Bst[j][k] for k in range(m))
            mu[i][j] = d / norm[j]
            for k in range(m): v[k] -= mu[i][j]*Bst[j][k]
        Bst.append(v)
        norm.append(sum(x*x for x in v))
    return Bst, mu, norm

def babai(cv, B):
    Bst, mu, norm = gso(B)
    cv = list(cv)
    n = len(B)
    # compute dot products with Bst
    for i in range(n-1, -1, -1):
        d = sum(cv[k]*Bst[i][k] for k in range(len(cv)))
        t = d / norm[i]
        ti = int(t.numerator // t.denominator + (1 if 2*abs(t.numerator % t.denominator) >= t.denominator else 0))
        # nearest integer
        ti = -(-t.numerator // t.denominator) if (t.numerator % t.denominator)*2 > t.denominator else t.numerator // t.denominator
        if ti:
            for k in range(len(cv)): cv[k] -= ti*Fraction(B[i][k])
    return cv

random.seed(3)
for it in range(6):
    cv = babai(cv, ZrF)
    random.shuffle(ZrF)
print("after babai: terms", nterms(cv), "cost", cost(cv))

# phase 2: +-1 elimination (integer t, denominators unchanged)
changed = True
rounds = 0
while changed and rounds < 200:
    changed = False; rounds += 1
    for zi, z in enumerate(Zr):
        for i in li5idx:
            if cv[i] and z[i] in (1, -1):
                t = cv[i] if z[i] == 1 else -cv[i]
                # apply cv -= t*z (zeroes cv[i] since z[i]==±1)
                sgn = 1 if z[i] == 1 else -1
                tt = cv[i]*sgn
                for k in range(len(cv)): cv[k] -= tt*Fraction(z[k])
                assert cv[i] == 0
                changed = True
print("after +-1 elim: terms", nterms(cv), "cost", cost(cv), "rounds", rounds)
cv = babai(cv, ZrF)
print("after babai2: terms", nterms(cv), "cost", cost(cv))

used = [(names[i], cv[i]) for i in range(1, len(cv)) if cv[i]]
dens = [x.denominator for _,x in used]; nums = [abs(x.numerator) for _,x in used]
print("BEST: terms:", len(used), " max|num| =", max(nums), " lcm den =", math.lcm(*dens), "=", sp.factorint(math.lcm(*dens)))
json.dump(dict(names=names, c=[str(x) for x in cv]), open(pipe.cand_dir(D)+"/rel_min_atom.json","w"))
for nm, cf in used: print(f"  {nm:>26s}  {cf}")
