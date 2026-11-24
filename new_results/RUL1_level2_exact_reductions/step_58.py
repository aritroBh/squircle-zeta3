import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
from fractions import Fraction as Fr
import rul1, rul1la, rul1_named, rul1_reduce
D = rul1_reduce.load()[3]
print("free:", D['free'])
print("fnames:", D['fnames'])
for nm, row in zip(D['fnames'], D['M']):
    print("M row", nm, {f: row.get(f, 0) for f in D['free']})
rref = [rul1la.Row(dict(d), dict(h), piv) for d, h, piv in D['rref']]
res, hist = rul1la.reduce_vector({(0,-1,-1): Fr(1)}, rref)
print("residue of (0,-1,-1):", res)
# check: residue should be a combo of free words with the same value
import mpmath as mp, rul1num
mp.mp.dps = 50
v_res = rul1num.combo_num(res, dps=45)
v_tgt = rul1num.g_word_over_0m1_num((0,-1,-1))
print("residue value", mp.nstr(v_res, 12), "target", mp.nstr(v_tgt, 12))
# print all w3 relations
for i,(fam,tag,r) in enumerate(D['rels']):
    print(i, fam, tag, r)
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=600)
print(r.stdout); print(r.stderr[-500:])
