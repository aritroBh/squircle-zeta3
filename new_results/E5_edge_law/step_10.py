import os, subprocess
d17='research/08_theory/D17_zeta5_polygons'
r=subprocess.run(['grep','-rl','960',d17],capture_output=True,text=True)
print(r.stdout)
r=subprocess.run(['grep','-rn','z2z3',d17],capture_output=True,text=True)
print(r.stdout[:4000])
