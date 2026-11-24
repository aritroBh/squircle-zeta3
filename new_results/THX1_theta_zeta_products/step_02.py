
import os
print(open('step_01.py').read())
base = '<project>/research'
for d in sorted(os.listdir(base)):
    print(d)
