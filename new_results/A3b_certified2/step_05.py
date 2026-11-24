
import os, flint
print("python-flint version:", flint.__version__)
p = '<project>/paper/squircle_note'
print("PAPER path exists:", os.path.exists(p))
print("qsqrt2 exists:", os.path.exists(os.path.join(p, 'scripts', 'qsqrt2.py')))
print("cwd:", os.getcwd())
