lead2 = '<project>/research/08_theory/LEAD2_leading_proof'
for f in ['step_01.py','step_02.py']:
    print('===', f)
    print(open(lead2+'/'+f).read())
print(os.path.exists(lead2+'/progress.log'))
import os
