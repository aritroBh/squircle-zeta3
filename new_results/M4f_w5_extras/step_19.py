WORK='<project>/research/08_theory/M4f_w5_extras'
import os
# The harness seems to pre-save steps; check what step_18 is
for f in ['step_06.py','step_18.py']:
    p=f'{WORK}/{f}'
    if os.path.exists(p):
        print('===',f, os.path.getsize(p)); print(open(p).read()[:300])
