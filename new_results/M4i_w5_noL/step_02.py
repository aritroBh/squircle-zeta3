import os
base = '<project>'
th = os.path.join(base, 'research/08_theory')
for d in sorted(os.listdir(th)):
    print(d)
print('---M4f---')
m4f = os.path.join(th, 'M4f_w5_extras')
for f in sorted(os.listdir(m4f)):
    p = os.path.join(m4f, f)
    print(f, os.path.getsize(p))
