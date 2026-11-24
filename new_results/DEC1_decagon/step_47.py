import subprocess
d='research/08_theory/OCT6_octagon_closed/'
out=subprocess.run(['sed','-n','1,120p',d+'step_19.py'],capture_output=True,text=True)
print(out.stdout)