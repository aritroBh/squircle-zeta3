import os
BASE = '<project>'
Q = os.path.join(BASE,'research/08_theory/A3c_claude_elementary/Q')
print(sorted(os.listdir(Q)))
src = open(os.path.join(Q,'cmzv.py')).read()
print(len(src))
print(src[:6000])
