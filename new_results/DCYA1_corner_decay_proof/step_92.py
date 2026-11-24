import subprocess
f='<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed/step_38.py'
out=subprocess.run(['cat',f],capture_output=True,text=True).stdout
print(out[6000:11000])
