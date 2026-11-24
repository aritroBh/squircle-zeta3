import os
ROOT = '<project>'
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
Q = os.path.join(ROOT,'research/08_theory/A3c_claude_elementary/Q')
src = open(os.path.join(Q,'cmzv.py')).read()
print(len(src))
print(src[:12000])
