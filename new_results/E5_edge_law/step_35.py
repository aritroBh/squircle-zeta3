import subprocess, os
# trivial sanity: does the venv python run a script cleanly?
open('/tmp/trivial.py','w').write("print('hello'); import sys; print(sys.path[:4])")
r=subprocess.run(['<project>/.venv/bin/python','/tmp/trivial.py'],capture_output=True,text=True,cwd='.')
print(r.stdout); print("ERR:", r.stderr[-500:])
# check for sitecustomize/usercustomize in the venv
import glob
print(glob.glob('<project>/.venv/lib/python*/site-packages/sitecustomize*'))
print(glob.glob('<project>/.venv/lib/python*/site-packages/usercustomize*'))
print(glob.glob('<project>/.venv/lib/python*/site-packages/*.pth'))
for p in glob.glob('<project>/.venv/lib/python*/site-packages/*.pth'):
    print(p, open(p).read()[:500])
