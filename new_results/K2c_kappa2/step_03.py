import os
base = '<project>/research/08_theory'
for d in ['K2c_kappa2','K2b_corner_phase','K2_corner_decay','D13d_decay_theorem']:
    p = os.path.join(base,d)
    print('===', d, '===')
    if os.path.isdir(p):
        for f in sorted(os.listdir(p)):
            print('  ', f)
    else:
        print('  not a dir:', os.path.exists(p))
