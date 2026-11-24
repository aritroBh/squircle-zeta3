
import os, shutil
BASE = '<project>'
WT2 = os.path.join(BASE, 'research/08_theory/WT2_weight2_all_k')
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')

# copy fastred.py locally, and read the truncated middle
shutil.copy(os.path.join(WT2,'fastred.py'), os.path.join(WORK,'fastred.py'))
src = open(os.path.join(WT2,'fastred.py')).read()
print(len(src))
# print the middle section that was truncated
print(src[3500:9500])
