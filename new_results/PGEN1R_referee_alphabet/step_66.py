import subprocess
# find PGEN1's actual regular-octagon delta computation (which coefficient table)
out = subprocess.run(["grep","-rn","L3d\|L3k\|3397","--include=step_*.py","../PGEN1_polygon_alphabet/step_44.py"],capture_output=True,text=True).stdout
print(out[:500])
src = open('../PGEN1_polygon_alphabet/step_44.py').read()
print(src[:4000])
