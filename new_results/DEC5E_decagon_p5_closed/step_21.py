import subprocess, os
base = '<project>'
out = subprocess.run(['grep','-rln','order5', base+'/paper'], capture_output=True, text=True)
print(out.stdout[:2000])
out = subprocess.run(['grep','-rln','73125', base+'/paper', base+'/research/08_theory'], capture_output=True, text=True)
print(out.stdout[:2000])