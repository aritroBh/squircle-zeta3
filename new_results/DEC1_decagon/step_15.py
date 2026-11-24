import os
d='research/08_theory/OCT6_octagon_closed/'
for f in sorted(os.listdir(d)):
    print(f)
print('---- oct5 subdir ----')
o5=os.path.join(d,'oct5')
if os.path.isdir(o5):
    for f in sorted(os.listdir(o5)): print(f)