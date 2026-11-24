import subprocess
# check import chains for precision traps
out = subprocess.run(['grep','-rn','import','<project>/research/08_theory/M4e_w5_lchi/w5_cobracket.py','<project>/research/08_theory/M4e_w5_lchi/w5_precision.py'],capture_output=True,text=True)
print(out.stdout)
out = subprocess.run(['grep','-rn','check_conversion\\|dilog_exact','<project>/research/08_theory/'],capture_output=True,text=True)
print(out.stdout[:2000])
