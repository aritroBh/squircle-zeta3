import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
src = open("step_02_varpi1.py").read()
src = src.replace('''from mpmath import mp, mpf, pi as mpi, sqrt as msqrt''',
'''from mpmath import mp, mpf, mpc, pi as mpi, sqrt as msqrt''')
src = src.replace('''print("c1        =", mp.nstr(mpf(str(sp.N(c1, 45))), 40))
print("varpi1    =", mp.nstr(mpf(str(sp.N(varpi1_paper, 45))), 40))''',
'''print("c1        =", mp.nstr(mpc(str(sp.N(c1, 45))), 40))
print("varpi1    =", mp.nstr(mpc(str(sp.N(varpi1_paper, 45))), 40))''')
open("step_02_varpi1.py","w").write(src)
r = subprocess.run(["<project>/.venv/bin/python","step_02_varpi1.py"],
                   capture_output=True, text=True, cwd=".", timeout=600)
print(r.stdout[-1200:], r.stderr[-1000:] if r.returncode else "")
