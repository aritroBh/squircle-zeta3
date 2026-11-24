cor = "<project>/research/08_theory/COR4_M3alpha_classical/"
step = '''"""step12_specializations.py -- Task 3: delta22(F) at special angles, corrected definition.
Per point: letter images + c0^e as exact field element; delta22 per c-exponent e is
rational; combined with c0^e in the field and tested for zero."""
import sys, json, re
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor4_sym import *

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
R3 = json.load(open(G1 + "hpl_fit_k3_(0,0,0).json"))
def parse(name, coef):
    m = re.match(r"c\\^(-?\\d+) H\\[(.*)\\]", name)
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
    num, den = (coef.split("/") + ["1"])[:2]
    return e, w, Fr(int(num), int(den))
w4 = [parse(n, c) for n, c in R3 if len(parse(n, c)[1]) == 4]

def F_by_e():
    out = {}
    for e, w, cf in w4:
        out.setdefault(e, []).append((w, cf))
    return out
FE = F_by_e()

def delta22_at(images, dim):
    """images: dict letter a in {0,1,-1} -> g-vector over local basis.
    Returns dict wedgepair -> dict e -> Fr (rational part for each c-exponent)."""
    out = {}
    for e, terms in FE.items():
        S = {}
        for w, cf in terms:
            word = tuple(images[a] for a in reversed(w))
            S[word] = S.get(word, 0) + cf
        for k, v in delta22L2(S, dim).items():
            d = out.setdefault(k, {})
            d[e] = d.get(e, 0) + v
    return {k: d for k, d in out.items() if any(d.values())}

results = {}

# --- c0 = 1/sqrt2, basis (2, 1+sqrt2) ---
h = Fr(1,2)
images = {0: (-h, Fr(0)), 1: (h, Fr(1)), -1: (-h, Fr(1))}
d = delta22_at(images, 2)
# c0^e in Q(sqrt2): e=-1: sqrt2 (0,1); e=0: 1; e=1: (0,1/2)
cpow = {-1: (Fr(0),Fr(1)), 0: (Fr(1),Fr(0)), 1: (Fr(0),h)}
tot = {}
for k, dd in d.items():
    acc = (Fr(0), Fr(0))
    for e, v in dd.items():
        c = cpow[e]
        acc = (acc[0] + v*c[0], acc[1] + v*c[1])
    if acc != (Fr(0),Fr(0)): tot[k] = acc
results["1/sqrt2 (alpha=pi/2)"] = tot
print("c = 1/sqrt2 (alpha = pi/2):  delta22 =", tot if tot else "0  [PREDICTED 0: paper Thm thm:W4R, M3 classical here]")

# --- c0 = 1/2, basis (2,3): [c]=(-1,0), [1-c]=(-1,0), [1+c]=(-1,1)
images = {0: (Fr(-1),Fr(0)), 1: (Fr(1),Fr(0)), -1: (Fr(-1),Fr(1))}
d = delta22_at(images, 2)
cpow = {-1: Fr(2), 0: Fr(1), 1: Fr(1,2)}
tot = {}
for k, dd in d.items():
    acc = sum((cpow[e]*v for e, v in dd.items()), Fr(0))
    if acc: tot[k] = acc
results["1/2 (alpha=2pi/3)"] = tot
print("c = 1/2     (alpha = 2pi/3): delta22 =", tot if tot else "0  [automatic: only 2 letters -> Lambda^2 B2 = 0]")

# --- c0 = sqrt3/2, basis (2,3,2+sqrt3): [c]=( -1, 1/2, 0), [1-c]=(-1,0,-1), [1+c]=(-1,0,1)
images = {0: (Fr(-1),h,Fr(0)), 1: (Fr(1),Fr(0),Fr(1)), -1: (Fr(-1),Fr(0),Fr(1))}
d = delta22_at(images, 3)
cpow = {-1: (Fr(0),Fr(2,3)), 0: (Fr(1),Fr(0)), 1: (Fr(0),Fr(1,2))}  # sqrt3/2 -> (0,1/2); 2/sqrt3 = 2sqrt3/3
tot = {}
for k, dd in d.items():
    acc = (Fr(0), Fr(0))
    for e, v in dd.items():
        c = cpow[e]
        acc = (acc[0] + v*c[0], acc[1] + v*c[1])
    if acc != (Fr(0),Fr(0)): tot[k] = acc
results["sqrt3/2 (alpha=pi/3)"] = tot
print("c = sqrt3/2 (alpha = pi/3):  delta22 =", tot if tot else "0")

# --- rational points ---
def rational_point(p, q, primes):
    # letters: [c]=[p/q], [1-c]=[(q-p)/q], [1+c]=[(p+q)/q] factored over primes
    import sympy as sp
    def fac(n):
        v = [Fr(0)]*len(primes)
        n = abs(n)
        for pr, ex in sp.factorint(n).items():
            v[primes.index(pr)] += ex
        return v
    def vec(num, den):
        a, b = fac(num), fac(den)
        return tuple(x - y for x, y in zip(a, b))
    images = {0: vec(p, q), 1: vec(q, q - p), -1: vec(p + q, q)}  # -[1-c] = [q/(q-p)]
    d = delta22_at(images, len(primes))
    c0 = Fr(p, q)
    cpow = {e: c0**e for e in (-1, 0, 1)}
    tot = {}
    for k, dd in d.items():
        acc = sum((cpow[e]*v for e, v in dd.items()), Fr(0))
        if acc: tot[k] = acc
    return tot

for p, q, primes, tag in [(1,3,[2,3],"1/3"), (2,5,[2,3,5,7],"2/5"), (3,7,[2,3,5,7],"3/7"),
                          (5,11,[2,3,11],"5/11"), (4,7,[2,3,7,11],"4/7"), (5,9,[2,3,7],"5/9")]:
    tot = rational_point(p, q, primes)
    results[f"{tag} (rational)"] = tot
    print(f"c = {tag:>5} (rational):      delta22 =", tot if tot else "0" + ("  [automatic: 2 letters]" if len(primes)==2 else ""))
'''
open(cor+"step12_specializations.py","w").write(step)
import subprocess
r = subprocess.run([sys.executable, cor+"step12_specializations.py"], capture_output=True, text=True, timeout=300)
print(r.stdout)
print(r.stderr[-800:] if r.returncode else "exit 0")
