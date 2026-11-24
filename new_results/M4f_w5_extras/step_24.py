import subprocess
r = subprocess.run(['<project>/.venv/bin/python','-c',
'''
import sys; sys.dont_write_bytecode=True
sys.path.insert(0,"<project>/research/08_theory/M4f_w5_extras")
import nice_search as NS
sol = NS.solve(list(NS.S0), "dbg")
kinds = {}
for name, s, rec in sol["basis"]:
    kinds.setdefault(rec[0], []).append(rec)
for k,v in kinds.items(): print(k, len(v), v[:2])
'''], capture_output=True, text=True)
print(r.stdout[-2000:], r.stderr[-1000:])
