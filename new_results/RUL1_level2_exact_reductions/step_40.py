import pathlib
p = pathlib.Path('rul1.py')
src = p.read_text()
old = """    cS = P.pop(S)
    assert cS == 1
    out = {}
    for T, cT in P.items():
        for U, cU in st_reg(T).items():
            out[U] = out.get(U, Fr(0)) - cT*cU
    return {k: v for k, v in out.items() if v}"""
new = """    cS = P.pop(S)
    assert cS in (1, 2)
    out = {}
    for T, cT in P.items():
        for U, cU in st_reg(T).items():
            out[U] = out.get(U, Fr(0)) - cT*cU/cS
    return {k: v for k, v in out.items() if v}"""
assert old in src
src = src.replace(old, new)
p.write_text(src)
print("st_reg multiplicity fixed")

import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
import mpmath as mp
import rul1, rul1num
mp.mp.dps = 60
rel = rul1.relation_DIST(((2,1),(3,1)))
print("DIST((2,1),(3,1)):")
for w, c in sorted(rel.items()): print("   ", w, c)
print("residual:", mp.nstr(abs(rul1num.combo_num(rel, dps=45)),3))
# hand-check: sum_delta zeta(2,d1;3,d2) = (1/8) zeta(2,3)
vals = {}
for d1 in (1,-1):
    for d2 in (1,-1):
        S = ((2,d1),(3,d2))
        vals[(d1,d2)] = rul1num.zeta_stuple_num(S, dps=45)
lhs = sum(vals.values())
rhs = mp.mpf(1)/8*rul1num.zeta_stuple_num(((2,1),(3,1)), dps=45)
print("lhs =", mp.nstr(lhs,15), " rhs =", mp.nstr(rhs,15), " diff =", mp.nstr(abs(lhs-rhs),3))
# E3 on the log-tail word of DUAL combo
for w in [(-1,-1,1), (-1,0,1), (0,-1,1)]:
    S = rul1.word_to_stuple(w)
    v = (-1)**len(S)*rul1num.zeta_stuple_num(S, dps=45)
    print("E3", w, S, mp.nstr(v,15))
print("st_reg((1,1),(1,1)) =", rul1.st_reg(((1,1),(1,1))))
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=900)
print(r.stdout); print(r.stderr[-800:])
