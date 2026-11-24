import os
base = '<project>'
# List 08_theory contents
for d in sorted(os.listdir(os.path.join(base, 'research/08_theory'))):
    print(d)
