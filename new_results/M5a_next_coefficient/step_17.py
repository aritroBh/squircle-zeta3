base = '<project>'
print(open(f'{base}/research/08_theory/M4a_decomposition/omega4.json').read())
# find which step created omega4.json
import subprocess
r = subprocess.run(['grep', '-l', 'omega4.json', 'research/08_theory/M4a_decomposition/step_*.py'], 
                   capture_output=True, text=True, cwd=base)
print(r.stdout)