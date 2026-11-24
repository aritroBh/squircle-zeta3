import subprocess
src='<project>/research/08_theory/WT2_weight2_all_k'
r = subprocess.run(['grep','-in','tower',f'{src}/progress.log'],capture_output=True,text=True)
print(r.stdout)
r = subprocess.run(['grep','-in','tail',f'{src}/progress.log'],capture_output=True,text=True)
print(r.stdout[-3000:])
