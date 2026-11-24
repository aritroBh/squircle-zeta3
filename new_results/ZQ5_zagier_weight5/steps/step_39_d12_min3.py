
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, math, pickle, random
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
N = len(cands)
c = [Fraction(-x, a0) for x in a]
Zf = [[Fraction(x) for x in z] for z in Z]
li5idx = list(range(2, 2+N))
Zi = [[int(x) for x in z] for z in Z]
R = flint.fmpz_mat(Zi).lll(gram='exact')
Zr = [[int(R[i,j]) for j in range(R.ncols())] for i in range(R.nrows()) if any(R[i,j] for j in range(R.ncols()))]
def cost(cv): return max(abs(x.numerator) for x in cv[1:])
def nterms(cv): return sum(1 for x in cv[1:] if x)

def capped_elim(cv, cap):
    cv = list(cv)
    changed = True; rounds = 0
    while changed and rounds < 300:
        changed = False; rounds += 1
        # find the +-1 elim that reduces terms and keeps cost <= cap, prefer biggest cost reduction
        bestmove = None
        for z in Zr:
            for i in li5idx:
                if cv[i] and z[i] in (1,-1):
                    sgn = 1 if z[i]==1 else -1
                    tt = cv[i]*sgn
                    co = max(abs((cv[k]-tt*z[k]).numerator) for k in range(1,len(cv)) if k!=i or True)
                    if co <= cap:
                        if bestmove is None or co < bestmove[0]:
                            bestmove = (co, z, i, tt)
        if bestmove:
            co, z, i, tt = bestmove
            for k in range(len(cv)): cv[k] -= tt*Fraction(z[k])
            assert cv[i] == 0
            changed = True
    return cv

frontier = []
for cap in (913, 1200, 2000, 5000, 20000, 200000):
    cv2 = capped_elim(c, cap)
    frontier.append((cap, nterms(cv2), cost(cv2)))
    print(f"cap={cap}: terms={nterms(cv2)} cost={cost(cv2)}")
# choose cap=5000 representative
cv = capped_elim(c, 5000)
used = [(names[i], cv[i]) for i in range(1, len(cv)) if cv[i]]
dens = [x.denominator for _,x in used]; nums = [abs(x.numerator) for _,x in used]
print("chosen: terms:", len(used), " max|num| =", max(nums), " lcm den =", math.lcm(*dens), "=", sp.factorint(math.lcm(*dens)))
json.dump(dict(names=names, c=[str(x) for x in cv]), open(pipe.cand_dir(D)+"/rel_min_atom.json","w"))
