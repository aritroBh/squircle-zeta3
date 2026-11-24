import subprocess, sys, datetime
code = r'''"""Task 3: convert gpl_w5/gpl_w4pp G-expressions to level-8 word vectors, exactly as
verify_D.rebuild_target does. Saves pickles; prints support sizes."""
import sys, pickle, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/research/08_theory/A3c_claude_elementary/Q")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from fractions import Fraction as F
import cmzv as cm
from cmzv import forms_path, dlog, atom, mul, add, POINTS
from qsqrt2 import QS2 as Q, A, B
import gexpr

def build(gex):
    lhs = [{}, {}]
    for (w, k), coefficient in gex.items():
        forms = [dlog(p) for p in w]
        term = add((forms_path(forms, POINTS[Q(1)], POINTS[cm.A]), 1),
                   (forms_path(forms, POINTS[Q(1)], POINTS[cm.B]), -1))
        for _ in range(k):
            term = mul(term, cm.atom(("log", Q(2))))
        for j, c in enumerate((coefficient.p, coefficient.q)):
            lhs[j] = add((lhs[j], 1), (term, c))
    return lhs

for name, fn in [("w4pp", gexpr.w4pp_g_expression), ("w5", gexpr.w5_g_expression)]:
    t0 = time.time()
    gex = fn()
    vec = build(gex)
    weights = {len(w) for d in vec for w in d}
    print(f"{name}: support (rat, sqrt2) = ({len(vec[0])}, {len(vec[1])}); "
          f"word weights present: {sorted(weights)}; ({time.time()-t0:.0f}s)", flush=True)
    assert weights <= ({4} if name == "w4pp" else {5}), weights
    with open(f"wv_{name}.pkl", "wb") as fh:
        pickle.dump([{w: str(c) for w, c in d.items()} for d in vec], fh)
    print(f"   saved wv_{name}.pkl")
'''
open("step_wordvectors.py", "w").write(code)
r = subprocess.run([sys.executable, "step_wordvectors.py"], capture_output=True, text=True, timeout=2700)
print(r.stdout)
print("ERR", r.stderr[-600:])
