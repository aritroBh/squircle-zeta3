import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
lines = open('av2core.py').read().split('\n')
for i in range(205, 230):
    print(i+1, lines[i])
