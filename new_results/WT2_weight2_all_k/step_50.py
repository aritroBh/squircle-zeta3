import sys, importlib, pickle, os
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import fastred; importlib.reload(fastred)
from fastred import *
import genintegrand, buildintegrand
importlib.reload(genintegrand); importlib.reload(buildintegrand)

record = pickle.load(open('record.pkl','rb'))

def gm_str(Gm):
    parts = []
    for r, c in sorted(Gm[0].items(), key=lambda rc: float(rc[0].sympy())):
        parts.append(f"[{c.sympy()}]@({r.sympy()})")
    for jj, e in enumerate(Gm[1], start=1):
        if not e.is_zero(): parts.append(f"[{e.sympy()}]/(v-i)^{jj}")
    if not Gm[3].is_zero(): parts.append(f"+const {Gm[3].sympy()}")
    return " ".join(parts) if parts else "0"

for key in [(3,0),(2,1),(4,0),(3,1),(2,2)]:
    print(f"===== Gm {key}")
    prev = None
    for k in range(sum(key), 11):
        rec = record.get((k, key))
        if rec is None or 'Gm' not in rec: continue
        s = gm_str(rec['Gm'])
        stable = "  (same as previous)" if s == prev else ""
        print(f"  k={k}: {s}{stable}")
        prev = s
