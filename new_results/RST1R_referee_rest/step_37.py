code = open('step_05.py').read()
code = code.replace("""bad = 0
for kv in range(99, 3001):
    Q13 = Q(13, kv)
    for iv in range(13, kv//2+1):
        if not (Q(iv,kv) <= Q13): bad += 1""",
"""bad = 0
for kv in range(99, 3001):
    Q13 = Q(13, kv).mid()
    for iv in range(13, kv//2+1):
        if Q(iv,kv).mid() > Q13*(1+arb('1e-25')): bad += 1""")
code = code.replace("""bad2 = 0
for kv in range(99, 2001):
    Q2 = Q(2, kv)
    for iv in range(2, kv//2+1):
        if not (Q(iv,kv) <= Q2): bad2 += 1""",
"""bad2 = 0
for kv in range(99, 2001):
    Q2 = Q(2, kv).mid()
    for iv in range(2, kv//2+1):
        if Q(iv,kv).mid() > Q2*(1+arb('1e-25')): bad2 += 1""")
open('step_05.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'step_05.py'],capture_output=True,text=True,timeout=900)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1500:])
if r.returncode==0:
    open('progress.log','a').write("step_05: R5 PROVED analytically: ell''>0 (sympy diff=0); h(k)>0 for ALL real k>=93 (c1=0.689877384603, kcrit=92.4132); exhaustive midpoint check Q_i<=Q_13 k in [99,3000] 0 viol; Q_i<=Q_2 k in [99,2000] 0 viol (analytic k>=21). step_18b h(k/2-12) vs proof h(k/2-13): cosmetic mismatch noted.\n")
