import sympy as sp
print([n for n in dir(sp) if 'subres' in n.lower()])
# fix script
src = open('run_subres.py').read()
src = src.replace("from sympy.polys.subresultants_qq_zz import subresultants\n","")
src = src.replace("subs = subresultants(C4_1, Q7_1, a5)","subs = sp.subresultants(C4_1, Q7_1, a5)")
open('run_subres.py','w').write(src)
import subprocess, time, os
PY = '<project>/.venv/bin/python'
logf = open('subres.log','w')
p = subprocess.Popen(['nice','-n','10',PY,'run_subres.py'], stdout=logf, stderr=subprocess.STDOUT,
                     cwd='<project>/research/08_theory/CHN3_p5_locus')
print("PID:", p.pid); time.sleep(30); print(open('subres.log').read())
