import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
from fractions import Fraction as Fr
import rul1_reduce, sympy as sp
D = rul1_reduce.load()[3]
Ms = sp.Matrix([[row.get(f, Fr(0)) for f in D['free']] for row in D['M']])
print("Ms =", Ms)
print("det =", Ms.det())
print("Ms.inv() =", Ms.inv())
print("stored Minv =", D['Minv'])
rv = [Fr(0), Fr(0), Fr(1,8)]
print("Ms.inv()*rv =", Ms.inv()*sp.Matrix(rv))
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=300)
print(r.stdout); print(r.stderr[-500:])
