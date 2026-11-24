import subprocess, os
# Find which script wrote Lchi8_5_2500.json and look at w6_verify.py + m5b_indep_check.py
r=subprocess.run(['grep','-rl','Lchi8_5_2500','<project>/research/08_theory/M5b_w6_closed_form/'],capture_output=True,text=True)
print(r.stdout)
print(open('w6_verify.py').read())
