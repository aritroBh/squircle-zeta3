import os
print(os.getcwd())
print(os.listdir('.'))
src='<project>/research/08_theory/WT2_weight2_all_k'
dst=os.getcwd()
import shutil
for f in ['fastred.py','genintegrand.py','buildintegrand.py','record2.pkl','binteg_15.pkl']:
    shutil.copy(os.path.join(src,f), os.path.join(dst,f))
print("copied:", os.listdir(dst))
