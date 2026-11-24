import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, math, time, random, pickle
import mpmath as mp
import sympy as sp
import w5lib as W
import pipe
D=12
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
LL = pickle.load(open(pipe.cand_dir(D)+"/atom_lll.pkl","rb"))
names = LL["names"]; a = LL["exact_rels"][0][1]; a0 = a[0]
Z = LL["exact_zeros"]
N = len(cands); M = len(names) - 2 - N
print("a0 =", a0, "=", sp.factorint(abs(a0)), " zeros:", len(Z), " N:", N, " M:", M)
c = [Fraction(-x, a0) for x in a]
assert c[0] == -1
Zf = [[Fraction(x) for x in z] for z in Z]
Zm = sp.Matrix(Zf)
rank = Zm.rank()
print("zero rank:", rank)
li5idx = list(range(2, 2+N))

def cost(cv):
    return max(abs(x.numerator) for x in cv[1:])

def solve_drop(S):
    Mm = sp.Matrix([[Zf[j][i] for j in range(len(Zf))] for i in S])
    if Mm.rank() < len(S): return None
    b = sp.Matrix([-c[i] for i in S])
    try: t = Mm.solve(b)
    except Exception: return None
    out = []
    for x in t:
        x = sp.Rational(x); out.append(Fraction(int(x.p), int(x.q)))
    return out

def apply_t(t):
    cv = list(c)
    for j, tj in enumerate(t):
        if tj:
            for i in range(len(cv)): cv[i] += tj*Zf[j][i]
    return cv

def nterms(cv): return sum(1 for x in cv[1:] if x)

best = None
t0=time.time()
random.seed(11)
for trial in range(400):
    S2 = random.sample(li5idx, rank)
    t = solve_drop(S2)
    if t is None: continue
    cv2 = apply_t(t)
    key = (nterms(cv2), cost(cv2))
    if best is None or key < best[0]:
        best = (key, S2, cv2)
        print(f"  trial {trial}: terms={key[0]} max|num|={key[1]}")
print(f"random done ({time.time()-t0:.0f}s)")
cv = best[2]
used = [(names[i], cv[i]) for i in range(1, len(cv)) if cv[i]]
dens = [x.denominator for _,x in used]; nums = [abs(x.numerator) for _,x in used]
print("BEST: terms:", len(used), " max|num| =", max(nums), " lcm den =", math.lcm(*dens), "=", sp.factorint(math.lcm(*dens)))
for nm, cf in used: print(f"  {nm:>26s}  {cf}")
json.dump(dict(names=names, c=[str(x) for x in cv], S=list(best[1])),
          open(pipe.cand_dir(D)+"/rel_min_atom.json","w"))
'''
open('steps/step_37_d12_min.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_37_d12_min.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_37_d12_min: D=12 minimisation\n")
