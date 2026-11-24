import subprocess
r = subprocess.run(['<project>/.venv/bin/python','-c',
'''
import sys, math; sys.dont_write_bytecode=True
sys.path.insert(0,"<project>/research/08_theory/M4f_w5_extras")
from flint import ctx
import nice_search as NS
from w5_search import evaluate
ctx.prec = 2000
sol = NS.solve(list(NS.S0), "dbg")
for name, s, rec in sol["basis"]:
    try:
        v = evaluate(rec, 100)
    except Exception as e:
        print("FAIL", name, rec, type(e).__name__, e)
        break
else:
    print("all ok")
'''], capture_output=True, text=True)
print(r.stdout[-1500:], r.stderr[-500:])
