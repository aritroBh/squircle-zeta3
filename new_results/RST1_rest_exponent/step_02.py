import os
os.chdir('<project>/research/08_theory')
for d in ['CX6_decay_rest','DCY2_decay_first_correction','DCYA1_corner_decay_proof','RST1_rest_exponent','D13_Mk_asymptotics']:
    print('='*60)
    print(d)
    if os.path.isdir(d):
        for f in sorted(os.listdir(d)):
            print('   ', f)
