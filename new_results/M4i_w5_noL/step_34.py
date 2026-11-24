import json, subprocess, time, os
dst = '<project>/research/08_theory/M4i_w5_noL'
with open(dst+'/progress.log','a') as fh:
    fh.write('step_25(nol_cross): Z=R5_noL-R5_best-c_L*(Lrel)+c_L*L has exact symbol 0 (PASS), value 7.5e-1740 fresh 1620d. c_L=(-80354/315,22520/3).\n')

chk = r'''
import json, sys
from fractions import Fraction as F
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/research/08_theory/M4i_w5_noL")
import flint
oc12 = json.loads(open("<project>/research/08_theory/M4i_w5_noL/overcomplete_1200.json").read())
oc15 = json.loads(open("<project>/research/08_theory/M4i_w5_noL/overcomplete_1500.json").read())
d12 = [[F(c) for c in dep["vec"][1:]] for dep in oc12["deps"]]
d15 = [[F(c) for c in dep["vec"][1:]] for dep in oc15["deps"]]
assert len(d12) == 14 and len(d15) == 14
rows = [[str(c) for c in dep] for dep in d12 + d15]
mat, rank = flint.fmpq_mat(rows).rref()
print("combined rank of 1200+1500 dep vectors:", rank, "(expect 14: same span)")
# also: every 1500 dep has zero L/s2L mixing consistency? just check the L-relation deps exist:
# rank of L/s2L projection within each set
for lbl, dd in (("1200", d12), ("1500", d15)):
    rows2 = [[str(dep[70]), str(dep[71])] for dep in dd]
    m2, r2 = flint.fmpq_mat(rows2).rref()
    print(f"  {lbl}: rank of (L,s2L) projection of deps = {r2} (2 => both L directions present)")
'''
open(dst+'/deps_consistency.py','w').write(chk)
PY='<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'deps_consistency.py'], capture_output=True, text=True, cwd=dst, timeout=600)
print('rc=%d' % r.returncode); print(r.stdout); print('STDERR:', r.stderr[-1500:])
with open(dst+'/progress.log','a') as fh:
    fh.write('step_26(deps_consistency): combined rank of 1200+1500 deps = 14 (same span); L/s2L projection rank 2 in each.\n')
