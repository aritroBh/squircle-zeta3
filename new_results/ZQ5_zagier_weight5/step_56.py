import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, math, time, itertools, random
import mpmath as mp
import sympy as sp
import w5lib as W
import pipe
D=5
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
vals = pipe.load_vals(D, 1150)
LL = __import__('pickle').load(open(pipe.cand_dir(D)+"/atom_lll.pkl","rb"))
names = LL["names"]; a = LL["exact_rels"][0][1]; a0 = a[0]
Z = LL["exact_zeros"]
N = len(cands); M = len(names) - 2 - N
print("a0 =", a0, " zeros:", len(Z), " N:", N, " M:", M)
c = [Fraction(-x, a0) for x in a]
assert c[0] == -1
Zf = [[Fraction(x) for x in z] for z in Z]
Zm = sp.Matrix(Zf)
print("zero rank:", Zm.rank())
li5idx = list(range(2, 2+N))

def cost(cv):
    return max(abs(x.numerator) for x in cv[1:])

def solve_drop(S):
    Mm = sp.Matrix([[Zf[j][i] for j in range(len(Zf))] for i in S])
    b = sp.Matrix([-c[i] for i in S])
    if Mm.rank() < len(S): return None
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

best = None
# strategy (iii): top-numerator 22
S0 = sorted(li5idx, key=lambda i: -abs(c[i].numerator))[:len(Zf)]
for S in [S0]:
    t = solve_drop(S)
    if t:
        cv = apply_t(t); best = (cost(cv), S, cv)
print("strategy top-num:", best[0] if best else None)
# greedy pivot
S = []
t0=time.time()
while len(S) < len(Zf):
    cand_best = None
    for i in li5idx:
        if i in S or c[i] == 0: continue
        Mm = sp.Matrix([[Zf[j][k] for j in range(len(Zf))] for k in S+[i]])
        if Mm.rank() <= len(S): continue
        t = solve_drop(S+[i])
        if t is None: continue
        co = cost(apply_t(t))
        if cand_best is None or co < cand_best[0]:
            cand_best = (co, i)
    if cand_best is None: break
    S.append(cand_best[1])
t = solve_drop(S); cv = apply_t(t)
if best is None or cost(cv) < best[0]: best = (cost(cv), S, cv)
print(f"greedy pivot: {best[0]} ({time.time()-t0:.0f}s), |S|={len(S)}")
# random tries
random.seed(7)
t0=time.time()
for _ in range(150):
    S2 = random.sample(li5idx, len(Zf))
    t = solve_drop(S2)
    if t is None: continue
    cv2 = apply_t(t); co = cost(cv2)
    if co < best[0]: best = (co, S2, cv2)
print(f"after random: {best[0]} ({time.time()-t0:.0f}s)")
cv = best[2]
used = [(names[i], cv[i]) for i in range(1, len(cv)) if cv[i]]
dens = [x.denominator for _,x in used]; nums = [abs(x.numerator) for _,x in used]
print("terms:", len(used), " max|num| =", max(nums), " lcm den =", math.lcm(*dens), "=", sp.factorint(math.lcm(*dens)))
li5used = [u for u in used if u[0].startswith("Li5")]
print("Li5 terms:", len(li5used))
json.dump(dict(names=names, c=[str(x) for x in cv], S=list(best[1])),
          open(pipe.cand_dir(D)+"/rel_min_atom.json","w"))
'''
open('steps/step_31_d5_min2.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_31_d5_min2.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_31_d5_min2: D=5 atom minimisation\n")
