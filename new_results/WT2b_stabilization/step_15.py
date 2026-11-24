src='<project>/research/08_theory/WT2_weight2_all_k'
for f in ['step_45.py','step_46.py','step_47.py']:
    c = open(f'{src}/{f}').read()
    print('='*20, f, len(c))
    print(c[:3500])
