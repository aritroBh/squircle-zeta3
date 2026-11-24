import os, subprocess
paper = "<project>/paper/squircle_note/scripts/"
print("paper scripts exist:", os.path.exists(paper))
for f in ["ibp_reduce.py","polylog_exact.py","qsqrt2.py","w2_integrals.py","kummer.py"]:
    print("  ",f, os.path.exists(paper+f))
work = "<project>/research/08_theory/TOW1R_review/"
step = r'''
import sys, os
sys.dont_write_bytecode=True
WTS="<project>/research/08_theory/WTS_weight_k_structure/theory/"
PAPER="<project>/paper/squircle_note/scripts/"
sys.path.insert(0,PAPER); sys.path.insert(0,WTS)
from fractions import Fraction as Fr
import fastred
from fastred import *
from buildintegrand import build_integrand_pf
from ibp_reduce import reduce_k
def pf_is_zero(P): return not P[0] and all(e.is_zero() for e in P[1]) and all(e.is_zero() for e in P[2]) and P[3].is_zero()
for k in [3,4]:
    slow = reduce_k(k, verbose=False)
    fast = fastred.reduce_k_fast(k, integrand_pf=build_integrand_pf(k))
    # compare ALL canonical keys, not just degree 2
    allkeys = sorted(set(slow.keys()) | set(fast.keys()))
    bad=0
    for key in allkeys:
        sPF = fastred.pf_from_sympy(slow.get(key,0))
        fPF = fast.get(key, pf_zero())
        if not pf_is_zero(pf_sub(sPF,fPF)): bad+=1; print("  mismatch",k,key)
    print(f"k={k}: slow ibp_reduce == fast reducer, ALL {len(allkeys)} keys, {bad} bad")
'''
with open(work+"step_23.py","w") as f: f.write(step)
env=dict(os.environ); env["PYTHONDONTWRITEBYTECODE"]="1"
r=subprocess.run(["nice","-n","19","<project>/.venv/bin/python","step_23.py"],
                 cwd=work,capture_output=True,text=True,env=env,timeout=1200)
print(r.stdout); print("STDERR:",r.stderr[:2000])
with open(work+"progress.log","a") as fh: fh.write("step_23: CONTROL fast reducer == paper ibp_reduce, ALL keys k=3,4\n")
