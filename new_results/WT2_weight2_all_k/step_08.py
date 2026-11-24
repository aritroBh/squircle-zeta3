import os
base = '<project>'
with open(os.path.join(base, 'paper/ieee/main.tex')) as f:
    lines = f.readlines()
print(''.join(lines[1005:1030]))
