import os
base='<project>'
for root,dirs,files in os.walk(os.path.join(base,'research/08_theory')):
    level=root.replace(os.path.join(base,'research/08_theory'),'').count(os.sep)
    if level<=1:
        print(root)
        for f in files:
            print('   ',f)
