code = open('step_09.py').read()
code = code.replace("from mpmath import mpf, mpi if False else mpf\n", "from mpmath import mpf\n")
open('step_09.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'step_09.py'],capture_output=True,text=True,timeout=3600)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1500:])
