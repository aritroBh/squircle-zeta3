import os, subprocess
base='<project>'
# Check for existing A3h2 dir
print(subprocess.run(['ls',base+'/research/08_theory/'],capture_output=True,text=True).stdout)
print(open(base+'/research/08_theory/A3h_weight4_target/Q/relations.py').read())
