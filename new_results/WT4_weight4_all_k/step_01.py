import os, subprocess
base='<project>'
# list the theory folder
for d in ['research/08_theory']:
    p=os.path.join(base,d)
    print('==',p)
    try:
        for x in sorted(os.listdir(p)):
            print('  ',x)
    except Exception as e:
        print('ERR',e)
