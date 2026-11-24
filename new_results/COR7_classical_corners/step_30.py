import os
base = '<project>/research/08_theory/COR5_corner_classicality/'
print(os.listdir(base+'planner'))
for f in os.listdir(base+'planner'):
    print('='*20, f, '='*20)
    print(open(base+'planner/'+f).read()[:4000])
