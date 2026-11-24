import pathlib
p = pathlib.Path('rul1.py')
src = p.read_text()
src = src.replace("out = _addinto(out, stuple_set_to_words(st_reg(plain) if not convergent_stuple(plain) else {plain: Fr(1)}), -Fr(2**(d-w)))",
                  "out = _addinto(out, stuple_set_to_words(st_reg(plain) if not convergent_stuple(plain) else {plain: Fr(1)}), -Fr(1, 2**(w-d)))")
old = """    out = {(): Fr(1)}
    for a in reversed(w):
        new = {}
        for suff, c in out.items():
            for l, sgn in letter[a]:
                nw = (l,) + suff
                new[nw] = new.get(nw, Fr(0)) + c*sgn
        out = new
    return {k: v for k, v in out.items() if v}"""
new = """    out = {(): Fr(1)}
    for a in reversed(w):
        new = {}
        for pref, c in out.items():
            for l, sgn in letter[a]:
                nw = pref + (l,)
                new[nw] = new.get(nw, Fr(0)) + c*sgn
        out = new
    return {k: v for k, v in out.items() if v}"""
assert old in src
src = src.replace(old, new)
p.write_text(src)
print("patched")

import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
import mpmath as mp
import rul1, rul1num
mp.mp.dps = 60
rel = rul1.relation_DIST(((2,1),(3,1)))
print("DIST((2,1),(3,1)) residual:", mp.nstr(abs(rul1num.combo_num(rel, dps=45)),3))
rel = rul1.relation_DUAL((0,1,1))
print("DUAL((0,1,1)) =", rel)
print("residual:", mp.nstr(abs(rul1num.combo_num(rel, dps=45)),3))
# all four families, one sample each at several weights
import random
random.seed(7)
for w in [2,3,4,5]:
    rels = rul1.relations_at_weight(w)
    byfam = {}
    for fam, tag, r in rels: byfam.setdefault(fam, []).append(r)
    line = "w=%d " % w
    for fam, rs in byfam.items():
        r = random.choice(rs)
        res = abs(rul1num.combo_num(r, dps=40))
        line += "%s(n=%d,res=%s) " % (fam, len(rs), mp.nstr(res,2))
    print(line)
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=900)
print(r.stdout); print(r.stderr[-800:])
