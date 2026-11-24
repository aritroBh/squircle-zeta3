
import os
for d in ['RUL1_level2_exact_reductions']:
    p = os.path.join('<project>/research/08_theory', d)
    for root, dirs, files in os.walk(p):
        for f in sorted(files):
            fp = os.path.join(root, f)
            print(fp, os.path.getsize(fp))
# K3
th = '<project>/research/08_theory'
for d in sorted(os.listdir(th)):
    if 'K3' in d:
        print("K3DIR:", d)
