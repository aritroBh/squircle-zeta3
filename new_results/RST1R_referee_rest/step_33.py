import re
code = open('step_02.py').read()
code = code.replace("""        for j in js:
            if j > 0:
                exact = -a*r + j*(rho/(rho+q/SQ2)).log()
            else:
                exact = -a*r
            Etil  = -a*r - SQ2*j*q/(2*r+1+q/SQ2)
            marg = Etil - exact
            if not (marg >= 0):
                bad += 1; print("VIOLATION", r.str(4), q.str(4), j, marg.str(6))
            worst_d = marg if worst_d is None else min(worst_d, marg)""",
"""        for j in js:
            if j == 0:
                continue   # j=0 is the equality case of (5)
            exact = -a*r + j*(rho/(rho+q/SQ2)).log()
            Etil  = -a*r - SQ2*j*q/(2*r+1+q/SQ2)
            marg = Etil - exact
            if not (marg > 0):
                bad += 1; print("VIOLATION", r.str(4), q.str(4), j, marg.str(6))
            worst_d = marg if worst_d is None else min(worst_d, marg)""")
open('step_02.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run([PY,'step_02.py'],capture_output=True,text=True)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1500:])
if r.returncode==0:
    open('progress.log','a').write("step_02: R1 verified: Pade g'=u^2/((1+u)(2+u)^2) (sympy); |s+q|>=rho+q/sqrt2 identity+stress OK; (5) margin>0 for j>=1 on 11x6x6 grid (j=0 equality). PASS\n")
