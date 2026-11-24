
import os, difflib
ROOT = 'research/08_theory/D13s_fix_recertify'
WORK = os.path.join(ROOT, 'D13s_work'); FIX = os.path.join(WORK, 'fixed')

for (orig, fixed) in [('D13g_all_N__step_38.py','step_38.py'), ('D13h_final_cell__step_10_J192.py','step_10_J192.py')]:
    a = open(os.path.join(WORK, orig)).read().splitlines()
    b = open(os.path.join(FIX, fixed)).read().splitlines()
    print("="*30, "DIFF", orig, "->", fixed, "="*30)
    for line in difflib.unified_diff(a, b, lineterm='', n=1):
        if line.startswith(('---','+++','@@')): print(line)
        elif line.startswith(('-','+')): print(line)
    print()
