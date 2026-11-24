# Step 7: independent control at k=5.  Compare the final pieces from the ORIGINAL sympy reducer
# (paper/squircle_note/scripts/ibp_reduce.py, copied to indep_k5/, output indep_k5/ibp_reduce_k5.out)
# against the fast exact reducer (final_5.pkl), key by key, by exact-rational-point evaluation.
import sys, os, re, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import fastred
from fastred import pf_to_sympy
vv = fastred.v
txt = open("indep_k5/ibp_reduce_k5.out").read().splitlines()
cur = None; orig = {}
for ln in txt:
    m = re.match(r"\s*weight (\d+)\s+l1\^(\d+) l2\^(\d+):", ln)
    if m:
        cur = (int(m.group(2)), int(m.group(3))); orig[cur] = 0; continue
    if cur is not None and ln.strip() and not ln.startswith("====="):
        orig[cur] += sp.sympify(ln.strip(), locals={"v": vv, "sqrt": sp.sqrt})
fin = pickle.load(open("final_5.pkl", "rb"))
print("keys: original", sorted(orig), " fast", sorted(fin))
assert set(orig) == set(fin), "key sets differ"
pts = [sp.Rational(3, 7), sp.Rational(5, 4), sp.Rational(17, 3), sp.Rational(-9, 5)]
worst = 0
for key in sorted(fin, key=lambda x: (sum(x), x)):
    d = orig[key] - pf_to_sympy(fin[key])
    mx = max(abs(sp.N(d.subs(vv, p), 60)) for p in pts)
    worst = max(worst, mx)
    print(f"  key {key} weight {sum(key)+1}: max |original - fast| over 4 test points = {sp.N(mx, 5)}")
print("WORST:", sp.N(worst, 5), " (60-digit evaluation of exact algebraic expressions)")
