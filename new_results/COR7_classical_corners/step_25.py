import subprocess, os
base = '<project>/research/08_theory/COR5_corner_classicality/'
# grep for basis definition in step files 60-78
for n in list(range(60,79)):
    f = base+f'step_{n}.py'
    if os.path.exists(f):
        src = open(f).read()
        if 'BASIS' in src or 'basis' in src:
            print('='*20, f'step_{n}.py', '='*20)
            print(src[:3000])
