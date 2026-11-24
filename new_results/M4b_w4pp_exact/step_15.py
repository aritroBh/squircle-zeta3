import os
ROOT = '<project>'
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
src = open(os.path.join(ROOT,'research/08_theory/A3_W4_proof/w4_symbol.py')).read()
i = src.find('def conjecture_terms')
print(src[i:i+3000])
