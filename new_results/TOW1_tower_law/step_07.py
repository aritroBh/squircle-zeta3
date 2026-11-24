import os
p = '../WTS_weight_k_structure/theory'
print(sorted(os.listdir(p)))
print(open(os.path.join(p,'result.md')).read()[:12000])
