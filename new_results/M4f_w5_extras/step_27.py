import subprocess
r = subprocess.run(['<project>/.venv/bin/python','-c',
'''
import sys; sys.dont_write_bytecode=True
sys.path.insert(0,"<project>/research/08_theory/M4f_w5_extras")
import mpmath as mp
import m4f_common as MC
import json
from fractions import Fraction as F
MC.set_prec(1200)
sol_names=None
import nice_search as NS
sol = NS.solve(list(NS.S0), "dbg")
names = sol["names"]
import pickle, pathlib
cache = pathlib.Path("<project>/research/08_theory/M4f_w5_extras/vals_cache_1200.pkl")
print("cache exists:", cache.exists())
if cache.exists():
    vals = pickle.loads(cache.read_bytes())
    K = json.loads(pathlib.Path("<project>/research/08_theory/M4f_w5_extras/nice_ksat_s0.json").read_text())["K"]
    kv = [mp.fsum(F(c)*v for c,v in zip(krow,vals)) for krow in K]
    for j,v in enumerate(kv):
        if abs(v) < mp.mpf(10)**(-1200):
            print(j, "ZERO", mp.nstr(abs(v),3))
    print("nonzero small ones:")
    for j,v in enumerate(kv):
        if abs(v) < 1 and abs(v) > 0:
            print(j, mp.nstr(abs(v),5))
'''], capture_output=True, text=True)
print(r.stdout, r.stderr[-800:])
