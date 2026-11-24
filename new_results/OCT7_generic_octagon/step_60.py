src = open("step_06.py").read()
# 1. unique monomials only
src = src.replace(
"""for (a,na),(b,nb) in itertools.product(logs, repeat=2): basis.append((a*b, na+nb))""",
"""for (a,na),(b,nb) in itertools.combinations_with_replacement(logs, 2): basis.append((a*b, na+nb))""")
src = src.replace(
"""for (a,na),(b,nb),(c,nc) in itertools.product(logs, repeat=3): basis.append((a*b*c, na+nb+nc))""",
"""for (a,na),(b,nb),(c,nc) in itertools.combinations_with_replacement(logs, 3): basis.append((a*b*c, na+nb+nc))""")
# 2. Fraction conversion in the 700-digit check
src = src.replace(
"""        tgt700 = Li21(*({'Lam1':(x1,x2/x1),'Lam2':(x2,x1/x2),'Lam3':(b1,b2/b1),'Lam4':(b2,b1/b2)}[nm]), )""",
"""        aa,cc = {'Lam1':(x1,x2/x1),'Lam2':(x2,x1/x2),'Lam3':(b1,b2/b1),'Lam4':(b2,b1/b2)}[nm]
        tgt700 = Li21(f2m(aa), f2m(cc))""")
open("step_06.py","w").write(src)
print("patched")
import subprocess, sys
r = subprocess.run([sys.executable, "step_06.py"], capture_output=True, text=True, cwd="research/08_theory/OCT7_generic_octagon", timeout=3600)
print(r.stdout)
if r.stderr: print("STDERR:", r.stderr[-800:])
