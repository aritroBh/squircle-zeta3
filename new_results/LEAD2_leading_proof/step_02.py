import os
lead1 = '<project>/research/08_theory/LEAD1_leading_genfun'
lead2 = '<project>/research/08_theory/LEAD2_leading_proof'
for p in [lead1, lead2]:
    print('===', p, os.path.exists(p))
    if os.path.exists(p):
        for f in sorted(os.listdir(p)):
            print('   ', f, os.path.getsize(os.path.join(p,f)))
