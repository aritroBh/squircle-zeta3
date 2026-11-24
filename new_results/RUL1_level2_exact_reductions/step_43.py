import pathlib
p = pathlib.Path('rul1.py')
src = p.read_text()
src = src.replace("assert cS in (1, 2)", "assert cS >= 1")
# dual_word: add (-1)^n factor in relation_DUAL (keep dual_word returning raw pullback)
old = '''def relation_DUAL(w):
    """duality for convergent w: G(w) = G(dual(w))."""
    out = {w: Fr(1)}
    for v, c in dual_word(w).items():
        out[v] = out.get(v, Fr(0)) - c
    return {k: x for k, x in out.items() if x}'''
new = '''def relation_DUAL(w):
    """duality for convergent w: G(w) = (-1)^{|w|} G(pullback(reversed w)) (t -> (1-t)/(1+t))."""
    sgn = (-1)**len(w)
    out = {w: Fr(1)}
    for v, c in dual_word(w).items():
        out[v] = out.get(v, Fr(0)) - sgn*c
    return {k: x for k, x in out.items() if x}'''
assert old in src
src = src.replace(old, new)
# also fix the docstring of dual_word
src = src.replace('"""t -> (1-t)/(1+t): letters pull back as 0 -> (1)-(-1), 1 -> (0)-(-1), -1 -> -(-1); reverse word."""',
                  '"""pullback letters of t -> (1-t)/(1+t), applied to reversed word; value sign (-1)^{|w|} handled in relation_DUAL."""')
p.write_text(src)
print("patched")

import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
import time
from fractions import Fraction as Fr
import rul1, rul1la, rul1num
import mpmath as mp
mp.mp.dps = 50
# DUAL spot checks
for w in [(-1,), (0,1), (0,1,1), (-1,-1,1)]:
    rel = rul1.relation_DUAL(w)
    res = abs(rul1num.combo_num(rel, dps=40)) if rel else 0
    print("DUAL", w, "residual", mp.nstr(res,3))
t0=time.time()
for w in range(1,5):
    rels = rul1.relations_at_weight(w)
    cols = sorted(rul1.conv_words_of_weight(w))
    rows = [rul1la.Row(r, {i: Fr(1)}) for i, (fam, tag, r) in enumerate(rels)]
    rref, piv = rul1la.eliminate(rows, cols)
    free = [c for c in cols if c not in piv]
    print("w=%d: rels %d, rank %d, words %d, codim %d (DG %s)  %.1fs" % (
        w, len(rels), len(piv), len(cols), len(free), [1,2,3,5][w-1], time.time()-t0), flush=True)
    print("   families:", {f: sum(1 for x in rels if x[0]==f) for f in ["DS","RDS","DIST","DUAL"]})
    print("   free:", free)
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=1500)
print(r.stdout); print(r.stderr[-800:])
