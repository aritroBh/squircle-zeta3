import subprocess, os
base = '<project>/research/08_theory'
k1d = f'{base}/K1d_engine_fixed'
print(sorted(os.listdir(k1d)))
r = subprocess.run(['grep','-l','Jtab', *[f'{k1d}/{f}' for f in sorted(os.listdir(k1d)) if f.endswith('.py')]], capture_output=True, text=True)
print("builders:", r.stdout)
