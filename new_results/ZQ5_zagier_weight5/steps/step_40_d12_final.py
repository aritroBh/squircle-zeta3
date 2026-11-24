
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
def cost(cv): return max(abs(x.numerator) for x in cv[1:])
def nterms(cv): return sum(1 for x in cv[1:] if x)

# regenerate the 88-term 913 representative (cap 913)
Zi = [[int(x) for x in z] for z in Z]
R = flint.fmpz_mat(Zi).lll(gram='exact')
Zr = [[int(R[i,j]) for j in range(R.ncols())] for i in range(R.nrows()) if any(R[i,j] for j in range(R.ncols()))]
def capped_elim(cv, cap):
    cv = list(cv); changed = True; rounds = 0
    while changed and rounds < 300:
        changed = False; rounds += 1
        bestmove = None
        for z in Zr:
            for i in li5idx:
                if cv[i] and z[i] in (1,-1):
                    sgn = 1 if z[i]==1 else -1
                    tt = cv[i]*sgn
                    co = max(abs((cv[k]-tt*z[k]).numerator) for k in range(1,len(cv)))
                    if co <= cap and (bestmove is None or co < bestmove[0]):
                        bestmove = (co, z, i, tt)
        if bestmove:
            co, z, i, tt = bestmove
            for k in range(len(cv)): cv[k] -= tt*Fraction(z[k])
            changed = True
    return cv
cv913 = capped_elim(c, 913)
print("913-rep: terms", nterms(cv913), "cost", cost(cv913))

# try random full eliminations ranked by cost with a 1e6 ceiling
def solve_drop(S):
    Mm = sp.Matrix([[Zf[j][i] for j in range(len(Zf))] for i in S])
    if Mm.rank() < len(S): return None
    b = sp.Matrix([-c[i] for i in S])
    try: t = Mm.solve(b)
    except Exception: return None
    return [Fraction(int(x.p), int(x.q)) for x in (sp.Rational(y) for y in t)]
def apply_t(t):
    cv = list(c)
    for j, tj in enumerate(t):
        if tj:
            for i in range(len(cv)): cv[i] += tj*Zf[j][i]
    return cv
best = (cost(cv913), nterms(cv913), cv913)
random.seed(23)
import time
t0=time.time(); wins=0
for _ in range(250):
    S2 = random.sample(li5idx, len(Zf))
    t = solve_drop(S2)
    if t is None: continue
    cv2 = apply_t(t); co = cost(cv2)
    if co <= 1000000 and co < best[0]:
        best = (co, nterms(cv2), cv2); wins += 1
print(f"random full-elim better than 913 found: {wins} (best {best[0]}, terms {best[1]}) ({time.time()-t0:.0f}s)")
cv = cv913 if best[0] >= 913 else best[2]
print("CHOSEN: terms", nterms(cv), "cost", cost(cv))
json.dump(dict(names=names, c=[str(x) for x in cv]), open(pipe.cand_dir(D)+"/rel_min_atom.json","w"))
