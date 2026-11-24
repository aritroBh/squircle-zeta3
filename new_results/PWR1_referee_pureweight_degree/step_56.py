import os
cands = ['../rhombus_formula.py', '../rhombus_check.py', '../rhombus_check2.py']
for c in cands:
    print(c, os.path.exists(c), os.path.getsize(c) if os.path.exists(c) else '')
print(open('../rhombus_formula.py').read()[:3500])
