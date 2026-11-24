for d in ['WT2_weight2_all_k','WT2b_stabilization','CX4_wt2b_review','CX5_wt2_writeup','WT3_weight3_all_k','WTS_weight_k_structure']:
    import os
    p = os.path.join('..', d)
    print('='*20, d)
    for root, dirs, files in os.walk(p):
        depth = root[len(p):].count(os.sep)
        if depth > 1: continue
        for f in sorted(files):
            fp = os.path.join(root,f)
            print(f'{os.path.getsize(fp):>9} {fp}')
