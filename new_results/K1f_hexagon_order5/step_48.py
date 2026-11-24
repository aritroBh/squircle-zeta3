import subprocess, os
base = '<project>/research/08_theory'
k1e = f'{base}/K1e_hexagon_closed'
r = subprocess.run(['grep','-l','Jtab', *[f'{k1e}/step_{i:02d}.py' for i in range(1,22)]], capture_output=True, text=True)
print(r.stdout)
