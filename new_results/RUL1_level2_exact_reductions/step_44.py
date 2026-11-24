import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
import rul1
for fam, tag, r in rul1.relations_at_weight(2):
    print(fam, tag, "->", r)
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=300)
print(r.stdout); print(r.stderr[-800:])
