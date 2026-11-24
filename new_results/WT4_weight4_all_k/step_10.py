import os
base='<project>/research/08_theory'
src=open(os.path.join(base,'WT2_weight2_all_k/fastred.py')).read()
print(src[9000:18000])
