
import os, subprocess, time, re
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
tc = work + "/tow4_copy"
PY = base + "/.venv/bin/python"

# (a) direct sanity of the lcm-numerator checker: import from prove_E3b78
sanity = r'''
import sys; sys.path.insert(0, ".")
import sympy as sp
from prove_E3b78 import check_zero
p = sp.Symbol('p')
print("check_zero([p+1]) [must be False]:", check_zero([p+1]))
print("check_zero([p/(p+1), -p/(p+1)]) [must be True]:", check_zero([p/(p+1), -p/(p+1)]))
print("check_zero([1/(2*p-1), 1/(2*p-3), sp.Rational(-1,3)]) [must be False]:",
      check_zero([1/(2*p-1), 1/(2*p-3), sp.Rational(-1,3)]))
'''
r = subprocess.run([PY,"-B","-c",sanity], capture_output=True, text=True, cwd=tc)
print(r.stdout, r.stderr)

# (b) perturbed feeds MUST fail
src = open(tc + "/prove_adj78_fast.py").read()
# perturbation 1: flip sign of the diagonal direct feed term
p1 = src.replace("terms.append(p*tau2(delta-dp)*Ediag(dp, s)/(2*d-dp-1))",
                 "terms.append(-p*tau2(delta-dp)*Ediag(dp, s)/(2*d-dp-1))")
assert p1 != src
open(tc + "/perturb1.py","w").write(p1.replace("prove_adj78_fast", "perturb1"))
# perturbation 2: tau1(n) -> tau1(n)+1 in the A-feed only
p2 = src.replace("terms.append((p+1)*tau1(delta-dp)*eR(dp, s, 2*p, p+1)*rA0/(2*d-dp-1))",
                 "terms.append((p+1)*(tau1(delta-dp)+1)*eR(dp, s, 2*p, p+1)*rA0/(2*d-dp-1))")
assert p2 != src
open(tc + "/perturb2.py","w").write(p2.replace("prove_adj78_fast", "perturb2"))
# perturbation 3: fold constant i**(-6-delta) -> i**(-4-delta)
p3 = src.replace("binom * i**(-6-delta)", "binom * i**(-4-delta)")
assert p3 != src
open(tc + "/perturb3.py","w").write(p3.replace("prove_adj78_fast", "perturb3"))
for name in ["perturb1.py","perturb2.py","perturb3.py"]:
    r = subprocess.run(["nice","-n","10",PY,"-B",name,"7","2"], capture_output=True, text=True, cwd=tc, timeout=600)
    out = (r.stdout.strip().splitlines() or [""])[-1]
    failed_as_expected = ("== 0 : False" in r.stdout) or (r.returncode != 0)
    print(f"{name}: rc={r.returncode} last='{out[:90]}' -> {'CORRECTLY FAILS' if failed_as_expected else '*** VACUOUS PASS BUG ***'}")
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} fast-method vacuity tests: check_zero sane; 3 perturbed feeds all correctly FAIL\n")
