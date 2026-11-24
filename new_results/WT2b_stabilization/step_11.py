src='<project>/research/08_theory/WT2_weight2_all_k'
import subprocess
print(subprocess.run(['tail','-c','6000',f'{src}/progress.log'],capture_output=True,text=True).stdout)
