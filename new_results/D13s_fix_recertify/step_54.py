
import os, shutil
ROOT = 'research/08_theory/D13s_fix_recertify'
WORK = os.path.join(ROOT, 'D13s_work')
# merge progress log to root
src = open(os.path.join(WORK,'progress.log')).read()
with open(os.path.join(ROOT,'progress.log'),'a') as f:
    f.write(src)
print(open(os.path.join(ROOT,'progress.log')).read())
