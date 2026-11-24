import shutil, os
src = '<project>/research/08_theory/WT2_weight2_all_k'
dst = '<project>/research/08_theory/WT2b_stabilization'
for f in ['fastred.py','genintegrand.py','buildintegrand.py']:
    shutil.copy(os.path.join(src,f), os.path.join(dst,f))
print(open(dst+'/genintegrand.py').read())
