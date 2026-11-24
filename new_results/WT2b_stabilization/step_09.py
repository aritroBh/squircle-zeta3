src='<project>/research/08_theory/WT2_weight2_all_k'
for f in ['step_55.py','step_56.py','step_57.py','step_58.py','step_59.py','step_60.py']:
    print('='*20, f)
    print(open(f'{src}/{f}').read()[:3000])
