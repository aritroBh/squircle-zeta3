import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from qsqrt2 import QS2 as Q

wp = pickle.load(open("wt5_primitive.pkl","rb"))
FORMS = pickle.load(open("closed_forms_lower.pkl","rb"))

def fmt_qs2(z):
    p, q = z.p, z.q
    return "%s,%s" % (p, q)

out = ["# auto-generated coefficient tables for verify_oct5d3.py (fractions p,q = p+q*sqrt2)"]
for name in ['Wt1','Wt2','Wt3','Wt4']:
    out.append(f"{name} = {{")
    for lab, (a, b) in FORMS[name].items():
        out.append(f"    {lab!r}: ({a!r}, {b!r}),")
    out.append("}")
out.append("WT5_BETA = {")
for lab, z in wp['beta'].items():
    out.append(f"    {lab!r}: ({str(z.p)!r}, {str(z.q)!r}),")
out.append("}")
out.append("WT5_LI32 = {")
for lab, z in wp['li32'].items():
    out.append(f"    {lab!r}: ({str(z.p)!r}, {str(z.q)!r}),")
out.append("}")
out.append("WT5_PLANNER = {")
for nm, z in wp['planner'].items():
    out.append(f"    {nm!r}: ({str(z.p)!r}, {str(z.q)!r}),")
out.append("}")
open("coeff_tables.py","w").write("\n".join(out)+"\n")
print("coeff_tables.py written, %d lines" % len(out))
# quick sanity: import it
import subprocess
r = subprocess.run([sys.executable, "-c", "import sys; sys.path.insert(0,'.'); import coeff_tables as ct; print(len(ct.Wt1), len(ct.Wt2), len(ct.Wt3), len(ct.Wt4), len(ct.WT5_BETA), len(ct.WT5_LI32), len(ct.WT5_PLANNER))"], capture_output=True, text=True)
print(r.stdout, r.stderr[-300:] if r.returncode else "")
