import os, subprocess
base = '<project>'
print(os.listdir(base))
print('---')
for d in os.listdir(base):
    if 'ZQ5' in d:
        print(d, '->', sorted(os.listdir(os.path.join(base,d)))[:40])
