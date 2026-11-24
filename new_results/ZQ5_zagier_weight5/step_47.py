import subprocess, sys
r = subprocess.run(['nice','-n','10',sys.executable,'../M4i_w5_noL/lrel_verify.py'],
                   capture_output=True,text=True,timeout=2400, cwd='../M4i_w5_noL')
print(r.stdout[-3000:]); print("STDERR:", r.stderr[-2000:])
with open('progress.log','a') as f: f.write("step_23_theirlrel: ran M4i lrel_verify.py\n")
