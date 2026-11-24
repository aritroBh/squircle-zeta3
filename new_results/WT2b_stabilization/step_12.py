src='<project>/research/08_theory/WT2_weight2_all_k'
import subprocess
print(subprocess.run(['grep','-n','tail\|tower\|stabil\|higher part\|Higher\|principal'][:0]+['grep','-in','tail\|tower',f'{src}/progress.log'],capture_output=True,text=True).stdout[-5000:])
