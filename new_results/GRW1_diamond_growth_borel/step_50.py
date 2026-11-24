import subprocess, os, time
wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
code = r'''
import flint
print(flint.__version__ if hasattr(flint,'__version__') else "n/a")
print([a for a in dir(flint) if 'calc' in a.lower() or 'integ' in a.lower()])
print([a for a in dir(flint.acb) if 'int' in a.lower() or 'calc' in a.lower()])
print([a for a in dir(flint.arb) if 'int' in a.lower() or 'calc' in a.lower()])
'''
r = subprocess.run(["<project>/.venv/bin/python","-c", code], capture_output=True, text=True, cwd=wd)
print(r.stdout); print(r.stderr[:500])
