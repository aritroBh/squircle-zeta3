import os
base='<project>'
for sub in ['research/08_theory/WT2_weight2_all_k','research/08_theory/M4a_decomposition','research/08_theory/M4b2_w4pp_nice','research/08_theory/KSPAN_kernel_span','research/08_theory/WT3_weight3_all_k','research/08_theory/WTS_weight_k_structure']:
    p=os.path.join(base,sub)
    print('==',sub)
    try:
        for x in sorted(os.listdir(p)):
            print('  ',x)
    except Exception as e:
        print('ERR',e)
