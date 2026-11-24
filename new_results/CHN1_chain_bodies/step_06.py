import os
PILOT='<project>/research/02_ideas/lenses/experimental_pilots'
for f in ['chain4_coeffs.txt','chain4_generic_1_1_1_2_d105.txt','chain4_generic_3_1_2_5.txt']:
    p=os.path.join(PILOT,f)
    print('='*10,f)
    print(open(p).read()[:800])
