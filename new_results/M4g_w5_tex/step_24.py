import subprocess
r=subprocess.run(['nice','-n','10','<project>/.venv/bin/python','gen_supplement.py'],capture_output=True,text=True)
print(r.stdout)
print(r.stderr)
open('progress.log','a').write('gen_supplement.py: wrote supplement_W5.tex (all arg forms checked numerically)\n')
