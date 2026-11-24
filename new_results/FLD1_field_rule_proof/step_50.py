import subprocess
G2='<project>/research/08_theory/G2_M4alpha/'
r = subprocess.run(['ls','-la',G2],capture_output=True,text=True)
print(r.stdout)
# check for result/fit files
r2 = subprocess.run(['bash','-c',f'cat {G2}*.out 2>/dev/null | tail -40; echo ===; cat {G2}result.md 2>/dev/null | head -60'],capture_output=True,text=True)
print(r2.stdout)