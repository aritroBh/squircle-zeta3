import subprocess
base='<project>'
w=base+'/research/08_theory/A3h2_w4_missing'
for i in range(1,16):
    print(f'--- step_{i:02d}.py ---')
    print(open(f'{w}/step_{i:02d}.py').read())
