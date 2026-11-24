import os, shutil, textwrap
base='<project>/research/08_theory'
work=os.path.join(base,'WT4_weight4_all_k')
os.makedirs(work, exist_ok=True)
# copy infrastructure (read-only sources)
for f in ['fastred.py','buildintegrand.py','genintegrand.py']:
    shutil.copy(os.path.join(base,'WT2_weight2_all_k',f), os.path.join(work,f))
for k in range(3,15):
    src=os.path.join(base,'WT2_weight2_all_k',f'binteg_{k}.pkl')
    if os.path.exists(src):
        shutil.copy(src, os.path.join(work,f'binteg_{k}.pkl'))
    else:
        print('missing', src)
print(sorted(os.listdir(work)))
open(os.path.join(work,'progress.log'),'a').write('step_setup: copied fastred/buildintegrand/genintegrand + binteg_3..14 from WT2_weight2_all_k\n')
