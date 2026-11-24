import subprocess
PY = '<project>/.venv/bin/python'
print(subprocess.run([PY,'-c',r'''
import mpmath as mp, inspect
print(inspect.getsource(mp.richardson))
'''], capture_output=True, text=True).stdout)
