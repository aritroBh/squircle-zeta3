import os
BASE = '<project>'
TH = os.path.join(BASE,'research/08_theory')
for d in sorted(os.listdir(TH)):
    print(d)
print('---A3_W4_proof---')
for f in sorted(os.listdir(os.path.join(TH,'A3_W4_proof'))):
    print(f, os.path.getsize(os.path.join(TH,'A3_W4_proof',f)))
