
import re
with open('step_03.py') as f:
    txt = f.read()
txt = txt.replace("VA = [np.array(v,float) for v in VA]", "VA = [np.array([float(x) for x in v]) for v in VA]")
txt = txt.replace("VB=[np.array(v,float) for v in VB]", "VB=[np.array([float(x) for x in v]) for v in VB]")
with open('step_03.py','w') as f:
    f.write(txt)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_03.py'],
                   capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-3000:] if r.returncode else "exit 0")
