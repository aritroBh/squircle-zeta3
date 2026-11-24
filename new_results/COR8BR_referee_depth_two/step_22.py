import subprocess
r = subprocess.run(['find', '<project>/research', '-name', 'F_symbol.json'], capture_output=True, text=True)
print(r.stdout)
r = subprocess.run(['ls', '<project>/research/08_theory/'], capture_output=True, text=True)
print(r.stdout)
