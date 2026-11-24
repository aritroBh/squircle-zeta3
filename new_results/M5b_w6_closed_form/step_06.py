import os
base='<project>/research/08_theory'
d=os.path.join(base,'M5b_w6_closed_form')
for f in sorted(os.listdir(d)):
    if f.endswith('.py'):
        print('='*20,f)
        print(open(os.path.join(d,f)).read())
