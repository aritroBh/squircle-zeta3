import os
base='<project>'
p=os.path.join(base,'research/08_theory/WT4_weight4_all_k')
print(os.listdir(p))
for f in ['WT4_weight4_all_k.brief']:
    fp=os.path.join(base,'research/08_theory',f)
    print('==',f)
    print(open(fp).read()[:4000])
