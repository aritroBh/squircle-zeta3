import subprocess
base='<project>/research/08_theory/HEXC1_hexagon_certificate'
print(subprocess.run(['grep','-l','make_relations']+ [f'{base}/step_{i:02d}.py' for i in range(1,102)],capture_output=True,text=True).stdout)
print(subprocess.run(['grep','-n','make_relations\|lower\|quotient\|rank', f'{base}/step_100.py'],capture_output=True,text=True).stdout)
