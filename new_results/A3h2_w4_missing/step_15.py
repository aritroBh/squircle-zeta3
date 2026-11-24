import subprocess
base='<project>'
print(subprocess.run(['ls','-la',base+'/research/08_theory/A3h2_w4_missing'],capture_output=True,text=True).stdout)
print(subprocess.run(['cat',base+'/research/08_theory/A3h2_w4_missing.brief'],capture_output=True,text=True).stdout[:3000])
