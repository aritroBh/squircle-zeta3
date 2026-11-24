import subprocess, os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','b3run2.py','[16384]','x'],
                   capture_output=True, text=True, timeout=1750)
print(r.stdout); print("STDERR:", r.stderr[-800:] if r.returncode!=0 else "(ok)")
with open('progress.log','a') as f:
    f.write("b3run2 [16384]: DONE\n")
