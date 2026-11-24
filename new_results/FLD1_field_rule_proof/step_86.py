import subprocess
D='<project>/research/08_theory/FLD1R_referee_field_rule/'
r=subprocess.run(['ls','-la',D],capture_output=True,text=True)
print(r.stdout[:4000])