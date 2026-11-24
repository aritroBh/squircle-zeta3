import sys, importlib, pickle, os
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import fastred; importlib.reload(fastred)
from fastred import *

record = pickle.load(open('record.pkl','rb'))
def pf_str(P, label=""):
    parts = []
    for r, c in sorted(P[0].items(), key=lambda rc: float(rc[0].sympy())):
        parts.append(f"[{c.sympy()}]@({r.sympy()})")
    for jj, e in enumerate(P[1], start=1):
        if not e.is_zero(): parts.append(f"[{e.sympy()}]/(v-i)^{jj}")
    if not P[3].is_zero(): parts.append(f"+const {P[3].sympy()}")
    return " ".join(parts) if parts else "0"

for key in [(2,0),(1,1)]:
    print(f"########## key {key}")
    for k in range(sum(key), 8):
        rec = record[(k,key)]
        print(f" k={k} W: {pf_str(rec['W'])}")
        print(f"      H: {pf_str(rec['H'])}")
