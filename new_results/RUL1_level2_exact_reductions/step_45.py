import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
import time
from fractions import Fraction as Fr
import rul1, rul1la
for fams in [("DS","DIST","DUAL"), ("DS","DUAL"), ("DS","DIST")]:
    t0=time.time()
    line = str(fams)+": "
    for w in range(1,6):
        rels = rul1.relations_at_weight(w, families=fams)
        cols = sorted(rul1.conv_words_of_weight(w))
        rows = [rul1la.Row(r, {i: Fr(1)}) for i, (fam, tag, r) in enumerate(rels)]
        rref, piv = rul1la.eliminate(rows, cols)
        free = [c for c in cols if c not in piv]
        line += " w%d: rank %d/%d codim %d;" % (w, len(piv), len(cols), len(free))
    print(line, " (%.1fs)" % (time.time()-t0), flush=True)
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=2400)
print(r.stdout); print(r.stderr[-800:])
