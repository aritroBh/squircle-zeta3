import os
print(os.getcwd())
print(os.listdir('.'))
base = '<project>'
with open(os.path.join(base, 'paper/ieee/main.tex')) as f:
    lines = f.readlines()
print(''.join(lines[830:875]))
