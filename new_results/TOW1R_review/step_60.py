import subprocess, os
work = "<project>/research/08_theory/TOW1R_review/"
step = r'''
import sys, os
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from myengine import *
bad=0
for p in range(2,25):
    # CX4 I3 relations: rho_{p+1,p-1} = -2p rho_pp/(p+1) ; 2p i rho_pp = (4p-3) rho_{p,p-1}
    if rho(p+1,p-1,0) != rho(p,p,0).scale(Fr(-2*p,p+1)): bad+=1; print("rel1 fail",p)
    if rho(p,p,0).scale(2*p)*IQ != rho(p,p-1,0).scale(4*p-3): bad+=1; print("rel2 fail",p)
    # rho_pp purely imaginary
    if not (rho(p,p,0).p0==0 and rho(p,p,0).p1==0): bad+=1; print("rho_pp not pure imag",p)
print("CX4 I3 intermediate relations + rho_pp pure-imaginary, p=2..24:", "0 bad" if bad==0 else f"{bad} bad")
# I4 relation: hhat_p = 2 i tau1 rho_pp = (1-i) rho_pp
bad=0
for p in range(1,25):
    if rho(p,p,0)*IQ*tau1(0).scale(2) != rho(p,p,0)*(QSI(1)-IQ): bad+=1; print("I4rel fail",p)
print("CX4 I4 relation 2i tau1 rho_pp = (1-i)rho_pp, p=1..24:", "0 bad" if bad==0 else f"{bad} bad")
'''
with open(work+"step_27.py","w") as f: f.write(step)
env=dict(os.environ); env["PYTHONDONTWRITEBYTECODE"]="1"
r=subprocess.run(["nice","-n","19","<project>/.venv/bin/python","step_27.py"],
                 cwd=work,capture_output=True,text=True,env=env,timeout=300)
print(r.stdout); print("STDERR:",r.stderr[:1500])
with open(work+"progress.log","a") as fh: fh.write("step_27: CX4 I3/I4 intermediate relations verified\n")
