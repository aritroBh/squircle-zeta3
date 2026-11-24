src='<project>/research/08_theory/WT2_weight2_all_k'
for f in ['step_40.py','step_41.py','step_42.py','step_43.py','step_44.py']:
    c = open(f'{src}/{f}').read()
    print('='*20, f, len(c))
    print(c[:2500])
