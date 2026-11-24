import os
base='<project>/research/08_theory/WTS_weight_k_structure/w4'
for f in ['identify_w4.md','result_w4.md']:
    print('='*30, f)
    try:
        print(open(os.path.join(base,f)).read()[:6000])
    except Exception as e:
        print('ERR', e)
