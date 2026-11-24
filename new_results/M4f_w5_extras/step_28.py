import subprocess
r = subprocess.run(['<project>/.venv/bin/python','-c',
'''
import sys; sys.dont_write_bytecode=True
sys.path.insert(0,"<project>/research/08_theory/M4f_w5_extras")
import mpmath as mp, json, pickle, pathlib
from fractions import Fraction as F
import m4f_common as MC
MC.set_prec(1200)
HERE = pathlib.Path("<project>/research/08_theory/M4f_w5_extras")
ks = json.loads((HERE/"nice_ksat_s0.json").read_text())
K, names = ks["K"], ks["names"]
vals = pickle.loads((HERE/"vals_cache_1200.pkl").read_bytes())
print("K rows:", len(K), "names:", len(names))
for j in [17,27,30]:
    sup = [(names[i], c) for i,c in enumerate(K[j]) if c]
    v = mp.fsum(F(c)*vv for c,vv in zip(K[j],vals))
    print("row", j, "value", mp.nstr(v,6))
    print("   support:", sup)
'''], capture_output=True, text=True)
print(r.stdout, r.stderr[-500:])
