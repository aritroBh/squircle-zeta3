import subprocess
target = '<project>/research/08_theory/COR8B_depth_two_two_terms'
r = subprocess.run(['grep', '-rn', 'hpl_fit_k3', target, '--include=*.py'], capture_output=True, text=True)
print(r.stdout[:3000])
