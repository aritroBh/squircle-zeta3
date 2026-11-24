import os
BASE = '<project>'
VD = os.path.join(BASE,'research/08_theory/A3h5_exact_certificate')
print(sorted(os.listdir(VD)))
src = open(os.path.join(VD,'verify_D.py')).read()
i = src.find('def rebuild_target')
print(src[i:i+3000])
