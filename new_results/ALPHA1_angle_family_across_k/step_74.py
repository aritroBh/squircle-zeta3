dbg = r'''
import sys
sys.path.insert(0, '.')
exec(open('step_30.py').read().split('Jc = to_hpl(Jexpr)')[0])
print("Jexpr terms:", len(Jexpr))
from collections import Counter
print("letters used:", Counter(a for (w,k) in Jexpr for a in w))
# try conv on each word, catch errors
fails = []
for (w, k), v in Jexpr.items():
    try:
        conv(w)
    except Exception as e:
        fails.append((w, str(e)))
print("conv failures:", len(fails))
for w, e in fails[:10]:
    print("   ", w, e)
# check a specific term: conv(('P','0','0')) etc
for w in [('P',), ('Q',), ('Q2',), ('P','0'), ('0','P'), ('1','Q2')]:
    try:
        print(w, "->", conv(w))
    except Exception as e:
        print(w, "ERR", e)
'''
open('dbg.py','w').write(dbg)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'dbg.py'], capture_output=True, text=True)
print(r.stdout[:4000])
print(r.stderr[-1500:] if r.returncode else "RC=0")
