import os, pickle, glob
base = '<project>'
theory = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')
for K in [9, 12, 18]:
    with open(os.path.join(theory, f'kern_{K}.pkl'),'rb') as f:
        KK = pickle.load(f)
    byd = {}
    for key in KK:
        byd.setdefault(sum(key), []).append(key)
    print('k=',K, {d: sorted(byd[d], key=lambda t:-t[0]) for d in sorted(byd)})
