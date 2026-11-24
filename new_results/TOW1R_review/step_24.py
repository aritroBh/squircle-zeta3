import subprocess, os
work = "<project>/research/08_theory/TOW1R_review/"
step = r'''
import sys, os
sys.dont_write_bytecode = True
WTS = "<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0, WTS)
from fractions import Fraction as Fr
from fastred import *
from buildintegrand import build_integrand_pf

def claimed_R20(k):
    # {0:(1-sqrt2)/4, a:1/4, -a:(2k-7)/4, -b:-(k-2)/2}
    d = {}
    d[ZERO]  = QSI(Fr(1,4), Fr(-1,4))          # (1-sqrt2)/4
    d[AQ]    = QSI(Fr(1,4))
    d[-AQ]   = QSI(Fr(2*k-7,4))
    d[-BQ]   = QSI(Fr(-(k-2),2))
    return d
def claimed_R11(k):
    d = {}
    d[ZERO]  = QSI(Fr(2*k-5,4), Fr(1,4))       # (2k-5+sqrt2)/4
    d[-AQ]   = QSI(Fr(-(2*k-5),2))
    return d

for k in [3,4,5]:
    ipf = build_integrand_pf(k)
    fin = reduce_k_fast(k, ipf)
    P20 = fin[(2,0)][0]; P11 = fin[(1,1)][0]
    ok20 = all(P20.get(r,ZERO)==c for r,c in claimed_R20(k).items()) and all(c.is_zero() or r in claimed_R20(k) for r,c in P20.items())
    ok11 = all(P11.get(r,ZERO)==c for r,c in claimed_R11(k).items()) and all(c.is_zero() or r in claimed_R11(k) for r,c in P11.items())
    print(f"k={k}: R20 match={ok20}  R11 match={ok11}")
    if not ok20:
        print("  got R20:", {r.sympy():c.sympy() for r,c in P20.items()})
        print("  exp    :", {r.sympy():c.sympy() for r,c in claimed_R20(k).items()})
    if not ok11:
        print("  got R11:", {r.sympy():c.sympy() for r,c in P11.items()})
        print("  exp    :", {r.sympy():c.sympy() for r,c in claimed_R11(k).items()})
    # also confirm no higher part remains at degree-2 keys
    print("   higher part (2,0):", fin[(2,0)][1], fin[(2,0)][2], " const:", fin[(2,0)][3])
'''
with open(work+"step_01.py","w") as f: f.write(step)
env = dict(os.environ); env["PYTHONDONTWRITEBYTECODE"]="1"
r = subprocess.run(["nice","-n","19","<project>/.venv/bin/python","step_01.py"],
                   cwd=work, capture_output=True, text=True, env=env, timeout=600)
print(r.stdout); print("STDERR:", r.stderr[:2000])
with open(work+"progress.log","a") as fh: fh.write("step_01: sanity kernels k=3,4,5 vs claimed formula\n")
