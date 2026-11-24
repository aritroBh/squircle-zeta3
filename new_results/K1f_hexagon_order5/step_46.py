import subprocess, os
base = '<project>/research/08_theory'
k1e = f'{base}/K1e_hexagon_closed'
r = subprocess.run(['grep','-l','lA0', *[f'{k1e}/{f}' for f in sorted(os.listdir(k1e)) if f.endswith('.py')]], capture_output=True, text=True)
print(r.stdout)
# look at the assembly step (step_08 mentioned) and a check step
print(open(f'{k1e}/step_08.py').read()[:3000])
