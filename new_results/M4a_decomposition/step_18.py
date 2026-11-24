import os
BASE = '<project>'
A3B2 = os.path.join(BASE,'research/08_theory/A3b_certified2')
print(sorted(os.listdir(A3B2)))
print('='*40)
rm = open(os.path.join(A3B2,'result.md')).read()
print(rm[:5000])
