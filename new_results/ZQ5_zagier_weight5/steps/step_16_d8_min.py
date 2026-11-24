
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json, math, itertools, time
import mpmath as mp
import sympy as sp
import flint
import w5lib as W
import pipe

D = 8
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
vals = pipe.load_vals(D, 1150)

R = pipe.phase_lll(D, C, cands, K, vals, P=830, use_atoms=True)
names = R["names"]; N = len(cands)
print("zero relations (a0=0):")
for mx,a,res in R["zeros"]:
    print("  max|a|=", mx, " resid=", mp.nstr(res,2))
    print("   ", [(names[i],a[i]) for i in range(len(a)) if a[i]])
print("top relations (a0, max|a|):")
for mx,a,res in R["rels"][:6]:
    print("  a0=", a[0], " max|a|=", mx)
# two-precision certificate: (R1-R2) in span of zeros
R1 = pipe.phase_lll(D, C, cands, K, vals, P=520, use_atoms=True)
r1 = R1["rels"][0][1]; r2 = R["rels"][0][1]
Zl = [a for _,a,_ in R["zeros"]]
d12 = [r1[i]-r2[i] for i in range(len(r1))]
Zm = sp.Matrix(Zl)
print("two-precision: (R520-R830) in Q-span of zero rels:", sp.Matrix(Zl+[d12]).rank() == Zm.rank())
print("exact rels:", len(R["exact_rels"]), " exact zeros:", len(R["exact_zeros"]))

# coefficient vector c on names (T coeff = -1): T = sum_{i>=1} c[i]*v_i
a = R["exact_rels"][0][1]; a0 = a[0]
c = [Fraction(-x, a0) for x in a]
assert c[0] == -1
Z = [ [Fraction(x) for x in zz] for _,zz,_ in R["exact_zeros"] ]

def cost(cv):
    return max(abs(x.numerator) for x in cv[1:])  # skip T

def solve_drop(S):
    # solve t: (c + sum t_j Z_j)[i] == 0 for i in S
    M = sp.Matrix([[Z[j][i] for j in range(len(Z))] for i in S])
    b = sp.Matrix([-c[i] for i in S])
    if M.rank() < len(S): return None
    try:
        t = M.solve(b)
    except Exception:
        return None
    out = []
    for x in t:
        x = sp.Rational(x)
        out.append(Fraction(int(x.p), int(x.q)))
    return out

def apply_t(t):
    cv = list(c)
    for j, tj in enumerate(t):
        if tj:
            for i in range(len(cv)):
                cv[i] += tj*Z[j][i]
    return cv

best = (cost(c), None, c)
li5idx = list(range(2, 2+N))
t0=time.time()
for sz in (5,4,3,2,1):
    for S in itertools.combinations(li5idx, sz):
        t = solve_drop(list(S))
        if t is None: continue
        cv = apply_t(t)
        co = cost(cv)
        if co < best[0]:
            best = (co, S, cv)
print("minimise-by-dropping done", f"{time.time()-t0:.0f}s", "best max|num|=", best[0], "drop S=", best[1])
cv = best[2]
# report
used = [(names[i], cv[i]) for i in range(1, len(cv)) if cv[i]]
dens = [x.denominator for _,x in used]
nums = [abs(x.numerator) for _,x in used]
print("terms:", len(used), " max|num|=", max(nums), " lcm den=", math.lcm(*dens))
for nm, cf in used:
    print(f"  {nm:>26s}  {cf}")
json.dump(dict(names=names, c=[str(x) for x in cv]), open(pipe.cand_dir(D)+"/rel_min.json","w"))
