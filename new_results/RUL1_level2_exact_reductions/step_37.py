import pathlib
p = pathlib.Path('rul1.py')
src = p.read_text()
old_ds = '''def relation_DS(u, v):
    """finite double shuffle for convergent u,v: u ш v - st(u)*st(v) transported to words."""
    lhs = shuffle(u, v)
    Su, Sv = word_to_stuple(u), word_to_stuple(v)
    rhs = stuple_set_to_words(stuffle(Su, Sv))
    out = dict(lhs)
    for w, c in rhs.items():
        out[w] = out.get(w, Fr(0)) - c
    return {k: x for k, x in out.items() if x}'''
new_ds = '''def relation_DS(u, v):
    """finite double shuffle for convergent u,v: G(u ш v) = G(u)G(v) = (-1)^{du+dv} zeta(Su)zeta(Sv);
    stuffle terms T transported with sign (-1)^{du+dv-dT} (zeta(T) = (-1)^{dT} G(w(T)))."""
    lhs = shuffle(u, v)
    Su, Sv = word_to_stuple(u), word_to_stuple(v)
    duv = len(Su) + len(Sv)
    out = dict(lhs)
    for T, c in stuffle(Su, Sv).items():
        sgn = Fr((-1)**(duv - len(T)))
        wT = stuple_to_word(T)
        cw = {wT: Fr(1)} if convergent_word(wT) else sh_reg(wT)
        for w, c2 in cw.items():
            out[w] = out.get(w, Fr(0)) - c*sgn*c2
    return {k: x for k, x in out.items() if x}'''
assert old_ds in src
src = src.replace(old_ds, new_ds)
old_rds = '''def relation_RDS(S):
    """regularized double shuffle for divergent stuple S: sh_reg(w(S)) - st_reg(S) transported."""
    w = stuple_to_word(S)
    lhs = sh_reg(w)
    rhs = stuple_set_to_words(st_reg(S))
    out = dict(lhs)
    for ww, c in rhs.items():
        out[ww] = out.get(ww, Fr(0)) - c
    return {k: x for k, x in out.items() if x}'''
new_rds = '''def relation_RDS(S):
    """regularized double shuffle (IKZ comparison) for divergent S:
    G^sh(w(S)) = sum_T c_T (-1)^{dS-dT} G(w(T))  where st_reg(S) = sum c_T T."""
    w = stuple_to_word(S)
    out = dict(sh_reg(w))
    dS = len(S)
    for T, c in st_reg(S).items():
        sgn = Fr((-1)**(dS - len(T)))
        wT = stuple_to_word(T)
        cw = {wT: Fr(1)} if convergent_word(wT) else sh_reg(wT)
        for ww, c2 in cw.items():
            out[ww] = out.get(ww, Fr(0)) - c*sgn*c2
    return {k: x for k, x in out.items() if x}'''
assert old_rds in src
src = src.replace(old_rds, new_rds)
p.write_text(src)
print("patched")

import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
import mpmath as mp
import rul1, rul1num
mp.mp.dps = 60
rel = rul1.relation_DS((0,-1), (0,-1))
print("DS((0,-1),(0,-1)) residual:", mp.nstr(abs(rul1num.combo_num(rel, dps=45)),3))
rel = rul1.relation_RDS(((1,1),(2,1)))
print("RDS((1,1),(2,1)) =", rel)
print("RDS residual:", mp.nstr(abs(rul1num.combo_num(rel, dps=45)),3))
rel = rul1.relation_DS((0,1), (0,1))
print("DS((0,1),(0,1)) =", rel)
print("residual:", mp.nstr(abs(rul1num.combo_num(rel, dps=45)),3))
rel = rul1.relation_DIST(((2,1),(3,1)))
print("DIST((2,1),(3,1)) residual:", mp.nstr(abs(rul1num.combo_num(rel, dps=45)),3))
rel = rul1.relation_DUAL((0,1,1))
print("DUAL((0,1,1)) =", rel)
print("residual:", mp.nstr(abs(rul1num.combo_num(rel, dps=45)),3))
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=900)
print(r.stdout); print(r.stderr[-800:])
