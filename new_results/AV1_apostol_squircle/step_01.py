import os, subprocess, sys
os.chdir('<project>/research/08_theory/AV1_apostol_squircle')
os.makedirs('steps', exist_ok=True)
print(os.getcwd())
print(open('progress.log').read() if os.path.exists('progress.log') else 'no log yet')
import mpmath
print('mpmath', mpmath.__version__, 'sympy ok, numpy/scipy ok')
