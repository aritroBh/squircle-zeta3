base='<project>/research/08_theory/WT2_weight2_all_k'
for n in [52,53,54,55,56,57]:
    print(f'===== step_{n:02d}.py =====')
    print(open(f'{base}/step_{n:02d}.py').read())
