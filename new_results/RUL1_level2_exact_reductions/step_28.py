import subprocess
PY = '<project>/.venv/bin/python'
print(subprocess.run([PY,'-c',r'''
import mpmath as mp, inspect
import mpmath.calculus.extrapolation as ex
src = inspect.getsource(ex.Extrapolation.richardson)
print(src[:3000])
'''], capture_output=True, text=True).stdout)
